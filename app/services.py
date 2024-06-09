import requests
import asyncio
import aiohttp
import networkx as nx
import osm2geojson
from geopy.distance import geodesic
import geopandas as gpd
from scipy.spatial import KDTree
import json
import ujson
from fastapi import FastAPI

app = FastAPI()


async def fetch(url, params):
    async with aiohttp.ClientSession() as session:
        async with session.get(url, params=params) as response:
            return await response.json()


async def get_roads_in_radius_osm(radius, start_point, end_point, center_point):
    overpass_url = "https://overpass-api.de/api/interpreter"

    async def fetch_roads(highway_type):
        overpass_query = f"""
            [out:json];
            (
                node(around: 1, {start_point[0]}, {start_point[1]});
                node(around: 1, {end_point[0]}, {end_point[1]});
                way[highway={highway_type}](around: {radius},{center_point[0]},{center_point[1]});
            );
            out geom;
        """
        return await fetch(overpass_url, params={"data": overpass_query})

    highway_types = [
        "primary",
        "primary_link",
        "secondary",
        "tertiary",
        "service",
        "trunk",
        "unclassified",
        "residental",
        "motorway"
    ]
    responses = await asyncio.gather(
        *[fetch_roads(highway_type) for highway_type in highway_types]
    )

    all_elements = []
    for response in responses:
        all_elements.extend(response["elements"])

    data = osm2geojson.json2geojson({"elements": all_elements})
    with open("main.json", "w") as f:
        json.dump(data, f, indent=2)


def get_charging_stations(radius, center_point):
    overpass_url = "https://overpass-api.de/api/interpreter"
    overpass_query = f"""
        [out:json];
        (
            node[amenity=charging_station](around: {radius},{center_point[0]},{center_point[1]});
        );
        out geom;
    """

    response = requests.get(overpass_url, params={"data": overpass_query})

    if response.status_code == 200:
        stations = []
        data = response.json()
        for station in data["elements"]:
            station_coords = (station["lat"], station["lon"])
            stations.append(station_coords)
        return stations
    else:
        raise Exception(f"Не удалось получить данные: {response.status_code}")


def find_nearest_road_point(road_graph, point):
    road_points = list(road_graph.nodes)
    kdtree = KDTree(road_points)
    distance, nearest_idx = kdtree.query(point)
    nearest_point = road_points[nearest_idx]
    return nearest_point


def build_graph(file_path, start_point, end_point, stations):
    data = gpd.read_file(file_path)
    road_graph = nx.DiGraph()

    for _, road in data.iterrows():
        road_geometry = road.geometry

        if (
            road_geometry.geom_type == "Point"
            or road_geometry.geom_type == "LineString"
        ):
            road_coordinates = list(road_geometry.coords)
        else:
            continue

        oneway = road['tags'].get('oneway', 'no')
        for i in range(len(road_coordinates) - 1):
            start = tuple(reversed(road_coordinates[i]))
            end = tuple(reversed(road_coordinates[i + 1]))
            distance = geodesic(start, end).meters
            if oneway == 'yes':
                road_graph.add_edge(start, end, length=distance)
            else:
                road_graph.add_edge(start, end, length=distance)
                road_graph.add_edge(end, start, length=distance)

    # Connect start_point to the nearest road point
    nearest_start = find_nearest_road_point(road_graph, start_point)
    road_graph.add_node(start_point)
    distance = geodesic(start_point, nearest_start).meters
    road_graph.add_edge(start_point, nearest_start, length=distance)
    road_graph.add_edge(nearest_start, start_point, length=distance)

    # Connect end_point to the nearest road point
    nearest_end = find_nearest_road_point(road_graph, end_point)
    road_graph.add_node(end_point)
    distance = geodesic(end_point, nearest_end).meters
    road_graph.add_edge(end_point, nearest_end, length=distance)
    road_graph.add_edge(nearest_end, end_point, length=distance)

    # Connect charging stations to the nearest road points
    for station in stations:
        nearest_station = find_nearest_road_point(road_graph, station)
        road_graph.add_node(station)
        distance = geodesic(station, nearest_station).meters
        road_graph.add_edge(station, nearest_station, length=distance)
        road_graph.add_edge(nearest_station, station, length=distance)

    return road_graph


def find_nearest_charging_station(point, stations):
    nearest_station = None
    min_distance = float("inf")

    for station in stations:
        distance = geodesic(point, station).meters
        if distance < min_distance:
            min_distance = distance
            nearest_station = station

    return nearest_station


# def optimize_route(graph, start, end, max_range, stations_on_range):
#     route = [start]
#     current_point = start
#     visited_stations = set()

#     while current_point != end:
#         try:
#             route_to_end = nx.astar_path(graph, source=current_point, target=end, weight="length")[1:]
#         except nx.NetworkXNoPath:
#             return None

#         if calc_distance_of_route(route_to_end) < max_range:
#             route += route_to_end
#             break
        
#         stations_on_range = [station for station in stations_on_range if station not in visited_stations]
#         if not stations_on_range:
#             return None
        
#         nearest_station = min(stations_on_range, key=lambda station: geodesic(end, station).meters)

#         if calc_distance_of_route([current_point, nearest_station]) > max_range:
#             nearest_station = min(stations_on_range, key=lambda station: geodesic(current_point, station).meters)
            
#         try:
#             route_segment = nx.astar_path(graph, source=current_point, target=nearest_station, weight="length")[1:]
#             route += route_segment
#             current_point = nearest_station
#             visited_stations.add(current_point)
#         except nx.NetworkXNoPath:
#             return None
        
#     return route

def optimize_route(graph, start, end, max_range, stations_on_range):
    route = [start]
    current_point = start
    visited_stations = set()
    old_range = max_range
    i = 0

    while current_point != end:
        try:
            route_to_end = nx.astar_path(graph, source=current_point, target=end, weight="length")[1:]
        except nx.NetworkXNoPath:
            print("Error with route")
            return None

        while max_range > 0:
            if route_to_end[i] != end:
                max_range -= calc_distance_of_route([route_to_end[i], route_to_end[i + 1]])
                i += 1
                route += [route_to_end[i]]
                current_point = route_to_end[i]
                nearest_station = min(stations_on_range, key=lambda station: calc_distance_of_route(nx.astar_path(graph, source=current_point, target=station, weight="length")))
                if max_range - calc_distance_of_route(nx.astar_path(graph, source=current_point, target=nearest_station, weight="length")) < (max_range / 10) or calc_distance_of_route(nx.astar_path(graph, source=current_point, target=nearest_station, weight="length")) < (max_range / 15) or (calc_distance_of_route(route_to_end) - calc_distance_of_route(route)) > max_range:
                    break
            else:
                current_point = route_to_end[i]
                break
        
        if current_point != end:
            stations_on_range = [station for station in stations_on_range if station not in visited_stations]
            if not stations_on_range:
                print("No stations")
                return None
                
            try:
                route_segment = nx.astar_path(graph, source=current_point, target=nearest_station, weight="length")[1:]
                route += route_segment
                current_point = nearest_station
                max_range = old_range
                i = 0
                visited_stations.add(current_point)
            except nx.NetworkXNoPath:
                print("Error with stations")
                return None
        
    return route


def calc_distance_of_route(route):
    result = 0
    for i in range(len(route) - 1):
        result += geodesic(route[i], route[i + 1]).meters
    return result


async def process_points(coords):
    start_point = (coords.start_lat, coords.start_lon)
    end_point = (coords.end_lat, coords.end_lon)
    center_point = (
        (start_point[0] + end_point[0]) / 2,
        (start_point[1] + end_point[1]) / 2)
    radius = geodesic(start_point, end_point).meters * 0.7

    await get_roads_in_radius_osm(radius, start_point, end_point, center_point)
    charging_stations = get_charging_stations(radius, center_point)
    print("Requst")
    road_graph = build_graph("main.json", start_point, end_point, charging_stations)
    print("Graph")

    optimal_route = optimize_route(road_graph, start_point, end_point, int(coords.max_distance), charging_stations)

    if optimal_route is None:
        print("Невозможно построить кратчайший путь.")

    print(calc_distance_of_route(optimal_route) / 1000, "km")
    return ujson.dumps({"route": optimal_route, 
                        "stations": charging_stations})
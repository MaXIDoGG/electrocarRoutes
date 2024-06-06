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
    road_graph = nx.Graph()

    for _, road in data.iterrows():
        road_geometry = road.geometry

        if (
            road_geometry.geom_type == "Point"
            or road_geometry.geom_type == "LineString"
        ):
            road_coordinates = list(road_geometry.coords)
        else:
            continue

        for i in range(len(road_coordinates) - 1):
            start = tuple(reversed(road_coordinates[i]))
            end = tuple(reversed(road_coordinates[i + 1]))
            distance = geodesic(start, end).meters
            road_graph.add_edge(start, end, length=distance)

    nearest_start = find_nearest_road_point(road_graph, start_point)
    road_graph.add_node(start_point)
    distance = geodesic(start_point, nearest_start).meters
    road_graph.add_edge(start_point, nearest_start, length=distance)

    nearest_end = find_nearest_road_point(road_graph, end_point)
    road_graph.add_node(end_point)
    distance = geodesic(start_point, nearest_end).meters
    road_graph.add_edge(end_point, nearest_end, length=distance)

    for station in stations:
        nearest_station = find_nearest_road_point(road_graph, station)
        road_graph.add_node(station)
        distance = geodesic(station, nearest_station).meters
        road_graph.add_edge(station, nearest_station, length=distance)

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


def optimize_route(graph, start, end, max_range):
    route = [start]
    current_point = start
    while current_point != end:
        route_to_end = nx.astar_path(
            graph, source=current_point, target=end, weight="length"
        )[1:]
        if calc_distance_of_route(route_to_end) < max_range:
            route += route_to_end
            break

        stations_on_range = get_charging_stations(max_range, current_point)
        if len(stations_on_range) == 0:
            return None
        need_station = min(stations_on_range, key=lambda st: geodesic(st, end).meters)

        route += nx.astar_path(
            graph, source=current_point, target=need_station, weight="length"
        )[1:]
        current_point = need_station

    return route


def calc_distance_of_route(route):
    result = 0
    for i in range(len(route) - 1):
        result += geodesic(route[i], route[i + 1]).meters
    return result


async def process_route(coords):
    start_point = (coords.start_lat, coords.start_lon)
    end_point = (coords.end_lat, coords.end_lon)
    center_point = (
        (start_point[0] + end_point[0]) / 2,
        (start_point[1] + end_point[1]) / 2,
    )
    radius = geodesic(start_point, end_point).meters * 0.7

    await get_roads_in_radius_osm(radius, start_point, end_point, center_point)
    charging_stations = get_charging_stations(radius, center_point)
    road_graph = build_graph("main.json", start_point, end_point, charging_stations)

    optimal_route = optimize_route(
        road_graph, start_point, end_point, coords.max_distance
    )

    if optimal_route is None:
        raise Exception("Невозможно построить кратчайший путь.")

    return ujson.dumps({"route": optimal_route})

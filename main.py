import requests
import networkx as nx
import osm2geojson
from geopy.distance import geodesic
import folium
import geopandas as gpd
from scipy.spatial import KDTree
import json
import tkinter as tk
from tkinter import messagebox

start_lat, start_lon = None, None
end_lat, end_lon = None, None
MAX_DISTANCE = None


def get_coordinates():
    global start_lat, start_lon
    global end_lat, end_lon
    global MAX_DISTANCE

    try:
        start_lat = float(entry_start_lat.get())
        start_lon = float(entry_start_lon.get())
        end_lat = float(entry_end_lat.get())
        end_lon = float(entry_end_lon.get())
        MAX_DISTANCE = float(entry_max_distance.get())

        # Проверка на максимальные и минимальные значения координат
        if (
            not (-90 <= start_lat <= 90)
            or not (-180 <= start_lon <= 180)
            or not (-90 <= end_lat <= 90)
            or not (-180 <= end_lon <= 180)
            or MAX_DISTANCE <= 0
        ):
            raise ValueError

        start_point = (start_lat, start_lon)
        end_point = (end_lat, end_lon)
        root.destroy()
        return start_point, end_point, MAX_DISTANCE
    except ValueError:
        messagebox.showerror("Ошибка", "Некорректные данные")
        return None, None, None


root = tk.Tk()
root.title("Введите координаты")

tk.Label(root, text="Начальная широта:").grid(row=0, column=0)
tk.Label(root, text="Начальная долгота:").grid(row=1, column=0)
tk.Label(root, text="Конечная широта:").grid(row=2, column=0)
tk.Label(root, text="Конечная долгота:").grid(row=3, column=0)
tk.Label(root, text="Максимальное расстояние (м):").grid(row=4, column=0)

entry_start_lat = tk.Entry(root)
entry_start_lat.grid(row=0, column=1)
entry_start_lon = tk.Entry(root)
entry_start_lon.grid(row=1, column=1)
entry_end_lat = tk.Entry(root)
entry_end_lat.grid(row=2, column=1)
entry_end_lon = tk.Entry(root)
entry_end_lon.grid(row=3, column=1)
entry_max_distance = tk.Entry(root)
entry_max_distance.grid(row=4, column=1)

tk.Button(root, text="OK", command=get_coordinates).grid(row=5, columnspan=2)

root.mainloop()

if (
    start_lat is None
    or start_lon is None
    or end_lat is None
    or end_lon is None
    or MAX_DISTANCE is None
):
    exit()

# Начальные и конечные координаты
START_POINT = (54.916267, 82.963470)
END_POINT = (55.038322, 82.942229)
CENTER_POINT = (
    (START_POINT[0] + END_POINT[0]) / 2,
    (START_POINT[1] + END_POINT[1]) / 2,
)


def get_roads_in_radius_osm(radius, start_point, end_point, center_point):
    overpass_url = "https://overpass-api.de/api/interpreter"
    overpass_query = f"""
        [out:json];
        (
            node(around: 1, {start_point[0]}, {start_point[1]});
            node(around: 1, {end_point[0]}, {end_point[1]});
            way[highway=primary](around: {radius},{center_point[0]},{center_point[1]});
            way[highway=primary_link](around: {radius},{center_point[0]},{center_point[1]});
            way[highway=secondary](around: {radius},{center_point[0]},{center_point[1]});
            way[highway=tertiary](around: {radius},{center_point[0]},{center_point[1]});
            way[highway=service](around: {radius},{center_point[0]},{center_point[1]});
            way[highway=trunk](around: {radius},{center_point[0]},{center_point[1]});
            way[highway=unclassified](around: {radius},{center_point[0]},{center_point[1]});
        );
        out geom;
    """
    response = requests.get(overpass_url, params={"data": overpass_query})

    if response.status_code == 200:
        print(200)
        data = response.json()
        data = osm2geojson.json2geojson(data)
        with open("main.json", "w") as f:
            json.dump(data, f, indent=2)
    else:
        messagebox.showerror(
            "Ошибка", f"Не удалось получить данные: {response.status_code}"
        )


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
        messagebox.showerror(
            "Ошибка", f"Не удалось получить данные: {response.status_code}"
        )
        return None


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

        stationsOnRange = get_charging_stations(max_range, current_point)
        if len(stationsOnRange) == 0:
            return None
        needStation = min(stationsOnRange, key=lambda st: geodesic(st, end).meters)

        route += nx.astar_path(
            graph, source=current_point, target=needStation, weight="length"
        )[1:]
        current_point = needStation

    return route


def calc_distance_of_route(route):
    result = 0
    for i in range(len(route) - 1):
        result += geodesic(route[i], route[i + 1]).meters
    return result


def visualize_route(route, start_point, end_point, center_point, stations):
    map_osm = folium.Map(location=center_point, zoom_start=12)

    folium.PolyLine(list(route), color="green", weight=5, opacity=1).add_to(map_osm)
    folium.CircleMarker(
        location=(start_point[0], start_point[1]), radius=7, color="blue", fill=False
    ).add_to(map_osm)
    folium.CircleMarker(
        location=(end_point[0], end_point[1]), radius=7, color="red", fill=False
    ).add_to(map_osm)
    for station in stations:
        folium.CircleMarker(
            location=station, radius=7, color="yellow", fill=False
        ).add_to(map_osm)

    map_osm.save("map_with_route.html")
    return map_osm


radius = geodesic(START_POINT, END_POINT).meters * 0.7  # В метрах
get_roads_in_radius_osm(radius, START_POINT, END_POINT, CENTER_POINT)
charging_stations = get_charging_stations(radius, CENTER_POINT)

road_graph = build_graph("main.json", START_POINT, END_POINT, charging_stations)
print(road_graph.nodes)
optimal_route = optimize_route(road_graph, START_POINT, END_POINT, MAX_DISTANCE)
if optimal_route == None:
    messagebox.showerror("Ошибка", "Невозможно построить кратчайший путь.")
    exit()
map_with_route = visualize_route(
    optimal_route, START_POINT, END_POINT, CENTER_POINT, charging_stations
)

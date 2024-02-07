import requests
import networkx as nx
from multiprocessing import Process
import osm2geojson
# import momepy
# from shapely.geometry import Polygon
from geopy.distance import geodesic
import folium
import geopandas as gpd
# import matplotlib.pyplot as plt
# import scipy as sp
# from scipy.spatial import KDTree
# import numpy as np
# import math
import json

# Начальные и конечные координаты (пример: Москва и Санкт-Петербург)
START_POINT = (55.031086, 82.921031)
END_POINT = (56.469176, 84.941186)
CENTER_POINT = ((START_POINT[0] + END_POINT[0]) / 2, (START_POINT[1] + END_POINT[1]) / 2)
MAX_DISTANCE = 150


# Запрос на OSM данных о дарогах в прямоугольнике
def get_roads_in_radius_osm(radius, start_point, end_point, center_point):
    overpass_url = "https://overpass-api.de/api/interpreter"
    overpass_query = f"""
        [out:json];
        (
            node(around: 100, {start_point[0]}, {start_point[1]});
            node(around: 100, {end_point[0]}, {end_point[1]});
            way[highway=primary](around: {radius},{center_point[0]},{center_point[1]});
            way[highway=secondary](around: {radius},{center_point[0]},{center_point[1]});
            
        );
        out geom;
    """

    # way[highway=tertiary](around: {radius},{center_point[0]},{center_point[1]});
    response = requests.get(overpass_url, params={'data': overpass_query})

    if response.status_code == 200:
        print(200)
        data = response.json()
        data = osm2geojson.json2geojson(data)
        with open("main.json", 'w') as f:
            json.dump(data, f, indent=2)
    else:
        return None


def get_charging_stations(radius, center_point):
    overpass_url = "https://overpass-api.de/api/interpreter"
    overpass_query = f"""
        [out:json];
        (
            node[amenity=charging_station](around: {radius},{center_point[0]},{center_point[1]});
        );
        out geom;
    """

    response = requests.get(overpass_url, params={'data': overpass_query})

    if response.status_code == 200:
        stations = []
        print(200)
        data = response.json()
        for station in data['elements']:
            station_coords = (station['lat'], station['lon'])
            stations.append(station_coords)
            print(station_coords)

        return stations
    else:
        return None


# Нахождение ближайшей точки дороги к заданной точке
def find_nearest_road_point(road_graph, point):
    nearest_point = None
    min_distance = float('inf')

    for road_point in road_graph.nodes:
        distance = geodesic(point, road_point).kilometers
        if distance < min_distance:
            min_distance = distance
            nearest_point = road_point

    return nearest_point


def build_graph(file_path, start_point, end_point, stations):
    # Загрузка данных из GeoJSON файла с помощью GeoPandas
    data = gpd.read_file(file_path)
    # Создание пустого графа
    road_graph = nx.Graph()

    # Обработка каждой дороги в данных
    for _, road in data.iterrows():
        # Получение геометрии дороги
        road_geometry = road.geometry

        # Получение координат точек, определяющих дорогу
        if road_geometry.geom_type == "Point" or road_geometry.geom_type == "LineString":
            road_coordinates = list(road_geometry.coords)
        else:
            continue

        # Добавление ребра в граф для каждой пары соседних точек
        for i in range(len(road_coordinates) - 1):
            start = tuple(reversed(road_coordinates[i]))
            end = tuple(reversed(road_coordinates[i + 1]))
            # Расчет дистанции (веса) между точками
            distance = geodesic(start, end).kilometers

            # Добавление ребра с весом в граф
            road_graph.add_edge(start, end, length=distance)

    nearest_start = find_nearest_road_point(road_graph, start_point)
    road_graph.add_node(start_point)
    distance = geodesic(start_point, nearest_start).kilometers
    road_graph.add_edge(start_point, nearest_start, length=distance)

    # Добавление конечной точки
    nearest_end = find_nearest_road_point(road_graph, end_point)
    road_graph.add_node(end_point)
    distance = geodesic(start_point, nearest_end).kilometers
    road_graph.add_edge(end_point, nearest_end, length=distance)

    for station in stations:
        nearest_station = find_nearest_road_point(road_graph, station)
        road_graph.add_node(station)
        distance = geodesic(station, nearest_station).kilometers
        road_graph.add_edge(station, nearest_station, length=distance)

    return road_graph


# Дистанция между вершинами графа
def calculate_distance(coord1, coord2):
    return geodesic(coord1, coord2).kilometers


def find_nearest_charging_station(point, stations):
    nearest_station = None
    min_distance = float('inf')

    for station in stations:
        distance = calculate_distance(point, station)
        if distance < min_distance:
            min_distance = distance
            nearest_station = station

    return nearest_station


# Нахождение кратчайшего пути
def optimize_route(graph, start, end, max_range, stations):
    start_station = find_nearest_charging_station(start, stations)
    end_station = find_nearest_charging_station(end, stations)
    route1 = nx.astar_path(graph, source=start, target=start_station, weight='length')
    route2 = nx.astar_path(graph, source=start_station, target=end_station, weight='length')[1:]
    route3 = nx.astar_path(graph, source=end_station, target=end, weight='length')[1:]
    route = route1 + route2 + route3
    return route


def visualize_route(route, start_point, end_point, center_point, stations):
    map_osm = folium.Map(location=center_point, zoom_start=12)

    folium.PolyLine(list(route), color="green", weight=5, opacity=1).add_to(map_osm)
    folium.CircleMarker(location=(start_point[0], start_point[1]), radius=7, color='blue', fill=False).add_to(map_osm)
    folium.CircleMarker(location=(end_point[0], end_point[1]), radius=7, color='red', fill=False).add_to(map_osm)
    for station in stations:
        folium.CircleMarker(location=station, radius=7, color='yellow', fill=False).add_to(
            map_osm)

    map_osm.save('map_with_route.html')
    return map_osm


# Радиус для запроса данных о дорогах
radius = calculate_distance(START_POINT, END_POINT) * 800  # В метрах
# Получение данных о дорогах
get_roads_in_radius_osm(radius, START_POINT, END_POINT, CENTER_POINT)
charging_stations = get_charging_stations(radius, CENTER_POINT)
# Построение графа
road_graph = build_graph("main.json", START_POINT, END_POINT, charging_stations)

# Применение алгоритма Дейкстры для поиска оптимального маршрута

optimal_route = optimize_route(road_graph, START_POINT, END_POINT, MAX_DISTANCE, charging_stations)

map_with_route = visualize_route(optimal_route, START_POINT, END_POINT, CENTER_POINT, charging_stations)

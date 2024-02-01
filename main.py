import requests
import networkx as nx
import osm2geojson
import momepy
from shapely.geometry import Polygon
from geopy.distance import geodesic
import folium
import geopandas as gpd
import matplotlib.pyplot as plt
import scipy as sp
from scipy.spatial import KDTree
import numpy as np
import math
import json

# Начальные и конечные координаты (пример: Москва и Санкт-Петербург)
start_point = (55.754410, 37.619114)
end_point = (55.557035, 37.708745)
center_point = ((start_point[0]+end_point[0])/2, (start_point[1]+end_point[1])/2)

# Запрос на OSM данных о дарогах в прямоугольнике
import requests
from pathlib import Path

def get_roads_in_radius_osm(radius, start_point, end_point, center_point):
    primary = "way[highway=primary](around:{},{},{});\n".format(radius, center_point[0], center_point[1])
    secondary = "way[highway=secondary](around:{},{},{});\n".format(radius, center_point[0], center_point[1])
    tertiary = "way[highway=tertiary](around:{},{},{});\n".format(radius, center_point[0], center_point[1])

    overpass_url = "https://overpass-api.de/api/interpreter"
    overpass_query = """
[out:json];
(
    node(around:10, {}, {});
    node(around:10, {}, {});
    {};
    {};
    {};
);
out geom;
    """.format(start_point[0], start_point[1], end_point[0], end_point[1], primary, secondary, tertiary)

    try:
        response = requests.get(overpass_url, params={'data': overpass_query})
        response.raise_for_status()
    except requests.exceptions.RequestException:
        return None
    
    data = response.json()
    with Path("main.json").open('w') as f:
        f.write(data)
    
def find_nearest_road_point(road_graph, point):
    nearest_point = min(road_graph.nodes, key=lambda x: geodesic(point, x).kilometers)
    return nearest_point

def build_graph(file_path, start_point, end_point):
    data = gpd.read_file(file_path)
    road_graph = nx.Graph()

    for _, road in data.iterrows():
        road_geometry = road.geometry
        if road_geometry.geom_type in ["Point", "LineString"]:
            road_coordinates = list(road_geometry.coords)
            for i in range(len(road_coordinates) - 1):
                start = tuple(reversed(road_coordinates[i]))
                end = tuple(reversed(road_coordinates[i + 1]))
                distance = geodesic(start, end).kilometers
                road_graph.add_edge(start, end, length=distance)

    nearest_start = find_nearest_road_point(road_graph, start_point)
    road_graph.add_node(start_point)
    distance = geodesic(start_point, nearest_start).kilometers
    road_graph.add_edge(start_point, nearest_start, length=distance)

    nearest_end = find_nearest_road_point(road_graph, end_point)
    road_graph.add_node(end_point)
    distance = geodesic(start_point, nearest_end).kilometers
    road_graph.add_edge(end_point, nearest_end, length=distance)

    return road_graph

# Дистанция между вершинами графа
def calculate_distance(coord1, coord2):
    return geodesic(coord1, coord2).meters


# Нахождение кратчайшего пути
def optimize_route(graph, start, end):
    return nx.shortest_path(graph, source=start, target=end, weight='length')

# Отображение маршрута на дороге
def visualize_route(route):
    map_osm = folium.Map(location=center_point, zoom_start=12)
    folium.PolyLine(list(road_graph), color="black", weight=5, opacity=1).add_to(map_osm)
    folium.PolyLine(list(route), color="green", weight=5, opacity=1).add_to(map_osm)
    folium.CircleMarker(location=(start_point[0], start_point[1]), radius=7, color='blue', fill=False).add_to(map_osm)
    folium.CircleMarker(location=(end_point[0], end_point[1]), radius=7, color='red', fill=False).add_to(map_osm)
    map_osm.save('map_with_route.html')
    return map_osm

def calc_all_distance(route):
    distance = 0
    for i in range(len(route) - 1):
        distance += calculate_distance(route[i], route[i + 1])
    return distance

# Радиус для запроса данных о дорогах
radius = calculate_distance(start_point, end_point) * 0.6 # В метрах
# Получение данных о дорогах
roads_data = get_roads_in_radius_osm(radius, start_point, end_point, center_point)

# Построение графа
road_graph = build_graph("main.json", start_point, end_point)
print("40%")
# Применение алгоритма Дейкстры для поиска оптимального маршрута

optimal_route = nx.shortest_path(road_graph, source=start_point, target=end_point, weight='length')
print("75%")

all_dis = calc_all_distance(optimal_route)
print(all_dis / 1000, "km")

map_with_route = visualize_route(optimal_route)
print("90%")

# Сохранение карты в файл (можно использовать для отображения в Jupyter Notebook)
map_with_route.save('map_with_route.html')
print("100%")
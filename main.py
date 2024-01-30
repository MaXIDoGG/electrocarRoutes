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
start_point = (55.031086, 82.921031)
end_point = (56.469176, 84.941186)
center_point = ((start_point[0]+end_point[0])/2, (start_point[1]+end_point[1])/2)

# Запрос на OSM данных о дарогах в прямоугольнике
def get_roads_in_radius_osm(radius, start_point, end_point):
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
    response = requests.get(overpass_url,params={'data': overpass_query})

    if response.status_code == 200:
        print(200)
        data = response.json()
        data = osm2geojson.json2geojson(data)
        with open("main.json", 'w') as f:
            json.dump(data, f, indent=2)
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

def build_graph(file_path, start_point, end_point):
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

    return road_graph

# Дистанция между вершинами графа
def calculate_distance(coord1, coord2):
    return geodesic(coord1, coord2).kilometers


# Нахождение кратчайшего пути
def optimize_route(graph, start, end):
    return nx.shortest_path(graph, source=start, target=end, weight='length')

def visualize_route(route):
    map_osm = folium.Map(location=center_point, zoom_start=12)
    
    # for i in range(len(route) - 1):
    #     node1 = route[i]
    #     node2 = route[i + 1]
    #     print(node1, " ", node2)
    folium.PolyLine(list(route), color="green", weight=5, opacity=1).add_to(map_osm)
        # i += 2
    folium.CircleMarker(location=(start_point[0], start_point[1]), radius=7, color='blue', fill=False).add_to(map_osm)
    folium.CircleMarker(location=(end_point[0], end_point[1]), radius=7, color='red', fill=False).add_to(map_osm)
    map_osm.save('map_with_route.html')
    return map_osm

# Радиус для запроса данных о дорогах
radius = calculate_distance(start_point, end_point) * 800 # В метрах
# Получение данных о дорогах
roads_data = get_roads_in_radius_osm(radius, start_point, end_point)
# Построение графа
road_graph = build_graph("main.json", start_point, end_point)
# print(road_graph.nodes)
# Применение алгоритма Дейкстры для поиска оптимального маршрута

optimal_route = nx.shortest_path(road_graph, source=start_point, target=end_point, weight='length')
# print(optimal_route)

map_with_route = visualize_route(optimal_route)

# Сохранение карты в файл (можно использовать для отображения в Jupyter Notebook)
# map_with_route.save('map_with_route.html')
from pydantic import BaseModel


class Coordinates(BaseModel):
    start_lat: float
    start_lon: float
    end_lat: float
    end_lon: float
    max_distance: float

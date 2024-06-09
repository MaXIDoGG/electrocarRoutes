from fastapi import APIRouter, HTTPException
from app.models import Coordinates
from app.services import process_points

router = APIRouter()


@router.post("/points")
async def get_points(coords: Coordinates):
    try:
        points_data = await process_points(coords)
        return points_data
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

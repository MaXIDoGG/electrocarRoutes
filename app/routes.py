from fastapi import APIRouter, HTTPException
from app.models import Coordinates
from app.services import process_route

router = APIRouter()


@router.post("/route")
async def get_route(coords: Coordinates):
    try:
        route_data = await process_route(coords)
        return route_data
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

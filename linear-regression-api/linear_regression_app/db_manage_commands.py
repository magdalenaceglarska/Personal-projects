from linear_regression_app import app, db
import json
from pathlib import Path
from linear_regression_app.models import PointsSet
from sqlalchemy.sql import text

@app.cli.group()
def db_manage():
    """Database management commands"""
    pass

@db_manage.command()
def add_data():
    """Add sample data to database"""
    try:
        pointssets_path = Path(__file__).parent / 'samples' / 'sets.json'
        with open(pointssets_path) as file:
            data_json = json.load(file)
        for item in data_json:
            pointsset = PointsSet(**item)
            db.session.add(pointsset)
        db.session.commit()
        print("Data has been added to database.")
    except Exception as exc:
        print("Unexpected error : {}.".format(exc))

@db_manage.command()
def remove_data():
    """Remove all data from database"""
    try:
        db.session.execute(text('TRUNCATE TABLE pointssets'))
        db.session.commit()
        print("Data has been removed from database.")
    except Exception as exc:
        print("Unexpected error : {}.".format(exc))
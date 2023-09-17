from flask import Flask
from config import Config
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy.sql import text
from flask_migrate import Migrate

app = Flask(__name__)
app.app_context().push()
app.config.from_object(Config)

db = SQLAlchemy(app)
migrate = Migrate(app, db)

#with app.app_context():
#    results = db.session.execute(text('show databases'))
#    for row in results:
#        print(row)

from linear_regression_app import points_sets
from linear_regression_app import models
from linear_regression_app import db_manage_commands
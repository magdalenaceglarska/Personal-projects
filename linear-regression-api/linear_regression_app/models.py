from linear_regression_app import db

class PointsSet(db.Model):
    __tablename__ = 'pointssets'
    id = db.Column(db.Integer, primary_key=True)
    x_values = db.Column(db.String(1000), nullable=False)
    y_values = db.Column(db.String(1000), nullable=False)

    def __repr__(self):
        return f'<{self.__class__.__name__}>: {self.x_values} {self.y_values}'
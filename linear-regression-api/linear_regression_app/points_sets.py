from linear_regression_app import app
from flask import jsonify

@app.route('/api/v1/points', methods=['GET'])
def get_points_sets():
    return jsonify({
        'success': True,
        'data': 'GET all points sets'
    })

@app.route('/api/v1/points/<int:set_id>', methods=['GET'])
def get_single_points_set(set_id: int):
    return jsonify({
        'success': True,
        'data': 'GET one points set'
    })

@app.route('/api/v1/points', methods=['POST'])
def create_points_set():
    return jsonify({
        'success': True,
        'data': 'New points set created'
    }), 201

@app.route('/api/v1/points/<int:set_id>', methods=['PUT'])
def update_points_set(set_id: int):
    return jsonify({
        'success': True,
        'data': f'Updated points set with id {set_id}'
    })

@app.route('/api/v1/points/<int:set_id>', methods=['DELETE'])
def delete_points_set(set_id: int):
    return jsonify({
        'success': True,
        'data': f'Deleted points set with id {set_id}'
    })
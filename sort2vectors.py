
def sort2vectors(v1,v2):

    points = zip(v1,v2)

    sorted_points = sorted(points)
    v1 = [point[0] for point in sorted_points]
    v2 = [point[1] for point in sorted_points]
    return v1,v2

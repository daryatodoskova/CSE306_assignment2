#define _CRT_SECURE_NO_WARNINGS 1
#include <iostream> 
#include <cmath>    
#include <random>
#include <string>
#include <algorithm>
#include <list>
#include <stdio.h>
#include <vector>
#include <limits>
#include <chrono>
#include "liblbfgs/lbfgs.c"
#include "polygon.cpp"
//#include "vector.cpp"
//#include "vector.h"

static std::default_random_engine engine(10); // random seed = 10
static std::uniform_real_distribution<double> uniform(0., 1.);

int main(){
    int n_vertices = 500;
    Polygon subjectPolygon = Polygon();
    Polygon clipPolygon = Polygon();
    Vector polygon1[4] = {Vector(0., 0.), Vector(0., 2.), Vector(2., 2.), Vector(2., 0.)};
    for (int j = 0; j < 4; j++)
    {
        clipPolygon.vertices.push_back(polygon1[j]);
    }

    Vector polygon2[3] = {Vector(2., 0.), Vector(1, 0.5), Vector(0.5, 0.7)};
   for (int j = 0; j < 3; j++)
    {
        subjectPolygon.vertices.push_back(polygon2[j]);
    }

    Polygon output = clipPolygonfunction(subjectPolygon, clipPolygon);
    std::vector<Polygon> polygons;
    polygons.push_back(output);
    polygons.push_back(clipPolygon);
    polygons.push_back(subjectPolygon);
    save_svg(polygons, "clip3.svg");
    
    for (int i = 0; i < n_vertices; i++){
        double x = uniform(engine);
        double y = uniform(engine);
        subjectPolygon.vertices.push_back(Vector(x, y));
    }
    
    double weights[subjectPolygon.vertices.size()];
    for (int i = 0; i < n_vertices; i++){
            weights[i] = 1;
    }

    std::vector<Polygon> result = voronoi(clipPolygon,subjectPolygon.vertices,weights);
    save_svg(result, "voronoi3.svg");

//----------------UNCOMMENT BELOW FOR POST-OPTIMISATION------------------

//     Vector polygon1[4] = {Vector(0., 0.), Vector(0., 1.), Vector(1., 1.), Vector(1., 0.)};
//     Polygon clipPolygon = Polygon();
//     for (int j = 0; j < 4; j++)
//     {
//         clipPolygon.vertices.push_back(polygon1[j]);
//     }

//     clipPolygon.area();
//     int n_vertices = 2000;

//     Polygon subjectPolygon = Polygon();
//     for (int i = 0; i < n_vertices; i++){
//         double x = uniform(engine);
//         double y = uniform(engine);
//         subjectPolygon.vertices.push_back(Vector(x, y));
//     }
    
//     double lambdas[subjectPolygon.vertices.size()];
//     double weights[subjectPolygon.vertices.size()];
//     Vector C = Vector(0.5,0.5);
//     double total;
//     for (int i = 0; i < n_vertices; i++){
//         Vector point = subjectPolygon.vertices[i];
//         Vector diff = C-point;
//         lambdas[i] = std::exp(-pow(diff.norm(), 2.) / 0.02);
//         total += lambdas[i];
//     }
//     for (int i = 0; i < n_vertices; i++)
//         {
//             lambdas[i] /= total;
//             weights[i] = 1;
//         }
//     OT optimizer = OT(clipPolygon,subjectPolygon.vertices, 100);
//     optimizer.execute();
//     std::vector<Polygon> result = optimizer.polygons;
//     result = voronoi(clipPolygon,subjectPolygon.vertices,weights);
    
//     result.push_back(subjectPolygon);
//     save_svg(result, "afterOptimisation5.svg");
}
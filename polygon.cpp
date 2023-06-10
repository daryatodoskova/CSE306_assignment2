#include "polygon.h"
#include <iostream>  
#include <cmath>     
#include <random>
#include <string>
#include <algorithm>
#include <list>
#include <stdio.h>
#include <vector>


// from https://pastebin.com/bEYVtqYy
// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
    void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
        FILE* f = fopen(filename.c_str(), "w+"); 
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        for (int i=0; i<polygons.size(); i++) {
            fprintf(f, "<g>\n");
            fprintf(f, "<polygon points = \""); 
            for (int j = 0; j < polygons[i].vertices.size(); j++) {
                fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
            }
            fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
            fprintf(f, "</g>\n");
        }
        fprintf(f, "</svg>\n");
        fclose(f);
    }


//----------------------------VECTOR------------------------------

Vector::Vector(double x, double y){
    coords[0] = x;
    coords[1] = y;
}

bool &Vector::operator==(const Vector &b){
    bool cond1 = coords[0] == b[0];
    bool cond2 = coords[1] == b[1];
    static bool res = cond1 && cond2;
    return res;
}

bool &Vector::operator!=(const Vector &b){
    bool cond1 = coords[0] != b[0];
    bool cond2 = coords[1] != b[1];
    static bool res = cond1 || cond2;
    return res;
}

Vector &Vector::operator+=(const Vector &b){
    coords[0] += b[0];
    coords[1] += b[1];
    return *this;
}

Vector &Vector::operator*=(const Vector &b){
    coords[0] *= b[0];
    coords[1] *= b[1];
    return *this;
}

Vector &Vector::operator/=(const Vector &b){
    coords[0] /= b[0];
    coords[1] /= b[1];
    return *this;
}

Vector &Vector::operator-=(const Vector &b)
{
    coords[0] -= b[0];
    coords[1] -= b[1];
    return *this;
}

const double &Vector::operator[](int i) const { return coords[i]; }
double &Vector::operator[](int i) { return coords[i]; }

Vector Vector::operator+(const Vector &a){
    return Vector(a[0] + coords[0], a[1] + coords[1]);
}

Vector Vector::operator+(const double a){
    return Vector(a + coords[0], a + coords[1]);
}

Vector Vector::operator-(const Vector &a){
    return Vector(coords[0] - a[0], coords[1] - a[1]);
}

Vector Vector::operator-(const double a){
    return Vector(coords[0] - a, coords[1] - a);
}

Vector Vector::operator*(const Vector &a){
    return Vector(a[0] * coords[0], a[1] * coords[1]);
}

Vector Vector::operator*(const double a){
    return Vector(a * coords[0], a * coords[1]);
}

Vector Vector::operator/(const Vector &a){
    return Vector(coords[0] / a[0], coords[1] / a[1]);
}

Vector Vector::operator/(const double a){
    return Vector(coords[0] / a, coords[1] / a);
}

double Vector::dot(const Vector &a){
    return a[0] * coords[0] + a[1] * coords[1];
}

double Vector::norm(){
    return sqrt(dot(*this));
}

Vector Vector::normalize(){
    return *this / norm();
}

//----------------------------EDGE------------------------------

Edge::Edge(const Vector &a, const Vector &b) {
    this->point_a = a;
    this->point_b = b;
};

//----------------------------POLYGON---------------------------

Polygon::Polygon(){
}

//compute polygon area
double Polygon::area()
{
    double res = 0;
    int n = vertices.size();
    if (n == 0)
    {
        return 0;
    }
    for (int i = 0; i < n; i++)
    {
        Vector p1 = vertices[i];
        Vector p2 = (i < n - 1) ? vertices[i + 1] : vertices[0];
        res += (p1[0] * p2[1] - p1[1] * p2[0]);
    }
    return std::abs(0.5 * res);
}

//---------------------------- UTILS -------------------------

//calculates the intersection point between an Edge and a line segment defined by two points a and b
Vector intersect(Vector &prevVertex, Vector &curVertex, Edge &clipEdge){
    Vector N = Vector(clipEdge.point_b[1] - clipEdge.point_a[1], clipEdge.point_a[0] - clipEdge.point_b[0]);
    double t = N.dot(clipEdge.point_a - prevVertex) / N.dot(curVertex - prevVertex);
    Vector P = prevVertex + (curVertex - prevVertex) * t;
    if (t < 0 || t > 1)
        return Vector(0., 0.);
    return P;
}

//checks if a given Vector point is inside a given Edge
bool inside(Vector &vertex, Edge &clipEdge)
{
    Vector N = Vector(clipEdge.point_b[1] - clipEdge.point_a[1], clipEdge.point_a[0] - clipEdge.point_b[0]) * (-1);
    bool test = N.dot(vertex - clipEdge.point_a) <= 0;
    return test;
}

//computes the area of a polygon given its vertices
double polygon_area(std::vector<Vector> vertices)
{
    double res = 0;
    int n = vertices.size();
    if (n == 0){return 0;}
    for (int i = 0; i < n; i++)
    {
        Vector p1 = vertices[i];
        Vector p2 = (i < n - 1) ? vertices[i + 1] : vertices[0];
        res += p1[0] * p2[1] - p1[1] * p2[0];
    }
    return std::abs(0.5 * res);
}

//---------------------------- Sutherland-Hodgman Algorithm -----------------------------------------
/* 
Clips the subjectPolygon by a convex clipPolygon (p.88 of lecture notes)
*/

Polygon clipPolygonfunction(Polygon &subjectPolygon, Polygon &clipPolygon)
{
    int prevIndex;
    for (int i = 0; i < clipPolygon.edges.size(); i++) //For each edge of the clip polygon
    {
        // Clip the subjectPolygon by a half-space
        Edge clipEdge = clipPolygon.edges[i];
        Polygon outPolygon = Polygon();
        for (int j = 0; j < subjectPolygon.vertices.size(); j++) // For each vertex of the subject polygon
        {
            // Test the subject polygon edge with vertices (j-1, j)
            Vector curVertex = subjectPolygon.vertices[j];
            if (j > 0)
                prevIndex = j - 1;
            else
                prevIndex = subjectPolygon.vertices.size() - 1;
            Vector prevVertex = subjectPolygon.vertices[prevIndex];
            // Compute inter between the infinite line supported by clipEdge and edge (j-1, j)
            Vector intersection = intersect(prevVertex, curVertex, clipEdge);
            if (inside(curVertex, clipEdge))
            {
                if (!inside(prevVertex, clipEdge))
                {
                    // The subject polygon edge crosses the clip edge, and we leave the clipping area
                    outPolygon.vertices.push_back(intersection);
                }
                outPolygon.vertices.push_back(curVertex);
            }
            else if (inside(prevVertex, clipEdge))
            {
                // The subject polygon edge crosses the clip edge, and we enter the clipping area
                outPolygon.vertices.push_back(intersection);
            }
        }

        subjectPolygon = outPolygon;
    }
    return subjectPolygon;
}



//---------------------------- VORONOI DIAGRAM --------------------------------

//clips a polygon by a line segment defined by two points a and b.
Polygon clipPolygonBis(Polygon &subjectPolygon, Vector M, Vector vectorIJ)
{
    int prevIndex;
    Polygon resPoly = Polygon();
    Vector vectorJI = vectorIJ * (-1);
#define inside(X) (X - M).dot(vectorIJ)
    for (int j = 0; j < subjectPolygon.vertices.size(); j++)
    {
        if (j<=0){
            prevIndex = subjectPolygon.vertices.size() - 1;
        }
        else{
            prevIndex = j - 1;
        }
        Vector curVertex = subjectPolygon.vertices[j];
        Vector prevVertex = subjectPolygon.vertices[prevIndex];
        double t = vectorJI.dot(M - prevVertex) / vectorJI.dot(curVertex - prevVertex);
        Vector intersection = (t >= 0 && t <= 1) ? prevVertex + (curVertex - prevVertex) * t : Vector(0., 0.);
        if (inside(curVertex) < 0)
        {
            if (!(inside(prevVertex) < 0))
            {
                resPoly.vertices.push_back(intersection);
            }
            resPoly.vertices.push_back(curVertex);
        }
        else if (inside(prevVertex) < 0)
        {
            resPoly.vertices.push_back(intersection);
        }
    }
    return resPoly;
}

//computes the Voronoi diagram given a clip polygon, a set of points, and their corresponding weights.
std::vector<Polygon> voronoi(Polygon &clipPolygon, std::vector<Vector> &points, const double *weights)
{
    Vector firstPoint, secondPoint, middlePoint, vectorIJ;
    double firstWeight, secondWeight;
    std::vector<Polygon> res;
    for (int i = 0; i < points.size(); i++)
    {
        firstWeight = weights[i];
        firstPoint = points[i];
        Polygon resPoly = clipPolygon;
        for (int j = 0; j < points.size(); j++)
        {
            secondPoint = points[j];
            secondWeight = weights[j];
            if (i == j){
                continue;
            }
            middlePoint = (firstPoint + secondPoint) * 0.5;
            vectorIJ = (secondPoint - firstPoint);
            middlePoint = middlePoint + vectorIJ * (firstWeight - secondWeight) / (2 * pow(vectorIJ.norm(), 2.));
            resPoly = clipPolygonBis(resPoly, middlePoint, vectorIJ);
        }
        res.push_back(resPoly);
    }
    return res;
}


//----------------------------  functions for OT -------------------------------------

//computes the triangulation of a polygon with respect to a given point
double triangulate(std::vector<Vector> vertices, Vector point)
{
    if (vertices.size() == 0){return 0;}
    Vector initial_vertex = vertices[0];
    double res = 0;
    for (int i = 0; i < vertices.size() - 2; i++)
    {
        Vector v1 = vertices[i + 1];
        Vector v2 = vertices[i + 2];
        std::vector<Vector> triangle;
        triangle.push_back(initial_vertex);
        triangle.push_back(v1);
        triangle.push_back(v2);
        double area = std::abs(polygon_area(triangle));
        res += (area / 6.) * ((initial_vertex - point).dot(initial_vertex - point) + (initial_vertex - point).dot(v1 - point) + (initial_vertex - point).dot(v2 - point) +
                              (v1 - point).dot(v1 - point) + (v1 - point).dot(v2 - point) +
                              (v2 - point).dot(v2 - point));
    }
    return std::abs(res);
}

//GD optimization to find optimal weights for the points in the given clip polygon
void gradient_descent(Polygon &clipPolygon, std::vector<Vector> &points, const double *lambdas, double eps, double step, double *weights)
{
    double fx, normGD;
    double gradient[points.size()], area[points.size()];
    std::vector<Polygon> polygons = voronoi(clipPolygon, points, lambdas);
    save_svg(polygons, "beforeOptimisation1.svg");
    double error = 1;
    int k = 1;
    while (error > eps)
    {
        error = 0;
        fx = 0;
        normGD = 0;
        polygons = voronoi(clipPolygon, points, weights);

        for (int i = 0; i < points.size(); i++)
        {
            area[i] = polygons[i].area();
            gradient[i] = (lambdas[i] - area[i]);
            normGD += gradient[i] * gradient[i];
        }
        for (int i = 0; i < points.size(); i++)
        {
            Vector point;
            std::vector<Vector> vertices;
            double temp;
            point = points[i];
            vertices = polygons[i].vertices;
            temp = (step / sqrt(normGD)) * gradient[i];
            error += temp * temp;
            weights[i] += temp;
            fx += (triangulate(vertices, point) - weights[i] * area[i] + lambdas[i] * weights[i]);
        }
        error = sqrt(error);
        k += 1;
        // if (k % 50 == 0){
        //     std::string filename = "file_" + std::to_string(k) + ".svg";
        //     save_svg(polygons, filename);
        //}
    }
}

//----------------------------OPTIMAL TRANSPORT--------------------------

OT::OT(Polygon &clipPolygon, std::vector<Vector> &points, int maximumiterations){
    this->clipPolygon = clipPolygon;
    this->points = points;
    this->max_iterations = maximumiterations;
}

OT::~OT(){
    if (m_x != NULL)
        {
            lbfgs_free(m_x);
            m_x = NULL;
        }
}

//executes the optimal transport calculation using the L-BFGS optimization algorithm.
int OT::execute(){
        lbfgsfloatval_t fx;
        int N = this->points.size();
        this->lambdas = (double *)malloc((N) * sizeof(double));

        lbfgsfloatval_t *m_x = lbfgs_malloc(N);
        if (m_x == NULL)
        {
            printf("ERROR: No memory allocated\n");
            return 1;
        }

        Vector C = Vector(0.5, 0.5);
        Vector diff;
        double tot = 0;
        double maximum = 0;
        for (int i = 0; i < N; i++)
        {
            lambdas[i] = std::exp(-pow(diff.norm(), 2.) / 0.02);
            tot += lambdas[i];
            diff = C - this->points[i];
        }
        for (int i = 0; i < N; i++)
        {
            m_x[i] = 1;
            lambdas[i] /= tot;
            if (lambdas[i] > maximum){
                maximum = lambdas[i];
            }
        }
        std::cout << "dirac:" << maximum << std::endl;
        this->polygons = voronoi(clipPolygon, this->points, lambdas);
        save_svg(polygons, "beforeOptimisation1.svg");
        lbfgs_parameter_t param;
        lbfgs_parameter_init(&param);
        param.max_iterations = this->max_iterations;

        int ret = lbfgs(N, m_x, &fx, _evaluate, _progress, this, &param);
        this->polygons = voronoi(clipPolygon, this->points, m_x);
        save_svg(polygons, "a.svg");
        printf("L-BFGS optimization terminated. Status code: %d\n", ret);
        free(lambdas);
        return ret;
}

//used by the L-BFGS optimizer to evaluate the function value and gradient at a given point
lbfgsfloatval_t OT::_evaluate(void *instance,const lbfgsfloatval_t *x,lbfgsfloatval_t *g,
                                    const int n, const lbfgsfloatval_t step){
        return reinterpret_cast<OT *>(instance)->evaluate(x, g, n, step);
}

//evaluates the function value and gradient at a given point during the L-BFGS optimization
lbfgsfloatval_t OT::evaluate(const lbfgsfloatval_t *x,lbfgsfloatval_t *g,const int n,const lbfgsfloatval_t step){
        lbfgsfloatval_t fx = 0.0;
        this->polygons = voronoi(clipPolygon, this->points, x);
        for (int i = 0; i < n; i++){
            std::vector<Vector> vertices = this->polygons[i].vertices;
            Vector point = this->points[i];
            double area = this->polygons[i].area();
            double temp = triangulate(vertices, point);
            fx += temp - x[i] * area + this->lambdas[i] * x[i];
            g[i] = area - this->lambdas[i];
        }
        fx = fx * -1.;
        return fx;
}

int OT::_progress(void *instance,const lbfgsfloatval_t *x,const lbfgsfloatval_t *g,const lbfgsfloatval_t fx,
                        const lbfgsfloatval_t xnorm,const lbfgsfloatval_t gnorm,const lbfgsfloatval_t step,
                        int n,int k,int ls){
        return reinterpret_cast<OT *>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
}

int OT::progress(const lbfgsfloatval_t *x,const lbfgsfloatval_t *g,const lbfgsfloatval_t fx,const lbfgsfloatval_t xnorm,
                const lbfgsfloatval_t gnorm,const lbfgsfloatval_t step,int n,int k,int ls){
        printf("Number of iteration is %d\n", k);
        printf(" Progress: fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
        printf("\n");
        return 0;   
}

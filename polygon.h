#define _CRT_SECURE_NO_WARNINGS 1
#include <iostream>
#include <vector>
#include "liblbfgs/lbfgs.h"

//----------------------------VECTOR-----------------------------

class Vector
{
public:
    explicit Vector(double x = 0., double y = 0.);
    bool &operator==(const Vector &b);
    bool &operator!=(const Vector &b);
    Vector &operator+=(const Vector &b);
    Vector &operator-=(const Vector &b);
    Vector &operator*=(const Vector &b);
    Vector &operator/=(const Vector &b);

    Vector operator+(const Vector &a);
    Vector operator+(const double a);

    Vector operator-(const Vector &a);
    Vector operator-(const double a);

    Vector operator*(const Vector &a);
    Vector operator*(const double a);

    Vector operator/(const Vector &a);
    Vector operator/(const double a);

    const double &operator[](int i) const;
    double &operator[](int i);

    double dot(const Vector &a);
    double norm();
    Vector normalize();

private:
    double coords[3];
};

//-----------------------------EDGE-------------------------------


class Edge{
    public:
        explicit Edge(const Vector &a,const Vector &b);
        Vector point_a;
        Vector point_b;
};

//--------------------------POLYGON-------------------------------

class Polygon {
    public:
        explicit Polygon();
        double area();
        double polygon_area(std::vector<Vector> vertices);
        std::vector<Edge> edges;
        std::vector<Vector> vertices;
        
};


//-------------------------OPTIMAL TRANSPORT----------------------

class OT{
public:
    std::vector<Vector> points;
    Polygon clipPolygon;
    double *weights;
    double *lambdas;
    int max_iterations;
    std::vector<Polygon> polygons;

    OT(Polygon &clipPolygon, std::vector<Vector> &points, int maxIterations);
    virtual ~OT();
    int execute();

protected:
    lbfgsfloatval_t *m_x;
    //constructs a power diagram whose weights are the variables passed in parameter to the "evaluate" function
    static lbfgsfloatval_t _evaluate(void *instance,const lbfgsfloatval_t *x,lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

    lbfgsfloatval_t evaluate(const lbfgsfloatval_t *x,lbfgsfloatval_t *g,const int n,const lbfgsfloatval_t step);

    static int _progress(void *instance,const lbfgsfloatval_t *x,const lbfgsfloatval_t *g,const lbfgsfloatval_t fx,
                        const lbfgsfloatval_t xnorm,const lbfgsfloatval_t gnorm,const lbfgsfloatval_t step,
                        int n,int k,int ls);

    int progress(const lbfgsfloatval_t *x,const lbfgsfloatval_t *g,const lbfgsfloatval_t fx,const lbfgsfloatval_t xnorm,
                const lbfgsfloatval_t gnorm,const lbfgsfloatval_t step,int n,int k,int ls);
    
};


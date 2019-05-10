#include <GL/glut.h>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <GL/glut.h>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>

#define ZERO 0.000000001

using namespace std;

class Vector;
class Point;
class Face;
class Mesh;

class Vector
{

public:
    double x, y, z;
    Vector()
    {
        x = y = z = 0;
    }
    Vector(double _x, double _y, double _z)
    {
        x = _x;
        y = _y;
        z = _z;
    }

    Vector(const Vector &vect)
    {
        x = vect.x;
        y = vect.y;
        z = vect.z;
    }

    Vector operator+(Vector &vect);
    Vector operator-(Vector &vect);
    Vector operator-() const;
    Vector operator*(double factor);
    Vector operator*(double factor) const;
    Vector operator/(double divisor);
    double operator*(Vector &vect);
    Vector operator|(Vector &vect);
    double len();
    Vector normalize();
    void print();
    bool isZero();
};

class Point : public Vector
{
public:
    Point() : Vector() {}

    Point(double _x, double _y, double _z) : Vector(_x, _y, _z){};

    Point(const Vector &&v)
    {
        x = v.x;
        y = v.y;
        z = v.z;
    }
};

class Face
{
public:
    vector<Point *> *mainVertice;
    Vector norm;
    vector<int> vertice;
    Mesh *mesh;

    float cs[4];
    float cd[4];
    float ca[4];
    float shininess = 100.0;

    Face()
    {
        mainVertice = NULL;
        mesh = NULL;
    }

    Face(vector<Point *> *mainVertices) : Face()
    {
        mainVertice = mainVertices;
    }

    Face(const Face &face)
    {
        mainVertice = face.mainVertice;
        norm = face.norm;
        vertice = face.vertice;
        mesh = face.mesh;
        for (int i = 0; i < 4; i++)
        {
            cs[i] = face.cs[i];
            cd[i] = face.cd[i];
            ca[i] = face.ca[i];
        }
        shininess = face.shininess;
    }
    void setMainVertices(vector<Point *> &mainVertices);
    void setAmbientColor(double red, double green, double blue, double alpha);
    void setSpecularColor(double red, double green, double blue, double alpha);
    void setDiffuseColor(double red, double green, double blue, double alpha);
    void setShininess(double value);
    void addPoint(int index);
    void recalNorm();
    void drawBorder();
    void drawColor();
    void setMesh(Mesh *m);
};

class Mesh
{
public:
    bool reverse = false;
    vector<Point *> mainVertice;
    vector<Face *> face;
    Mesh *dpO = NULL, *dpOX = NULL, *dpOY = NULL, *dpOZ = NULL;
    vector<Mesh *> dp, dpX, dpY, dpZ;
    Vector lightDirection;
    double csr, csg, csb;
    Mesh()
    {
        lightDirection = Vector(0, 0, 0);
        csr = 1, csg = 1, csb = 1;
    }
    Mesh(string filename) : Mesh()
    {
        // Read data from a file

        fstream f;
        f.open(filename);
        int nVertex = 0;
        int nFace = 0;
        f >> nVertex >> nFace;
        mainVertice.clear();
        for (int i = 0; i < nVertex; i++)
        {
            double _x, _y, _z;
            f >> _x >> _y >> _z;
            mainVertice.push_back(new Point(_x, _y, _z));
        }
        for (int i = 0; i < nFace; i++)
        {
            newFace();
            Face &lface = getLastFace();
            int vCount;
            f >> vCount;
            while (vCount--)
            {
                int x;
                f >> x;
                lface.addPoint(x);
            }
            double r, g, b, a;
            f >> r >> g >> b >> a;
            lface.setAmbientColor(r, g, b, a);
        }
    }
    void addPoint(Point &point);
    void addFace(Face &_face);
    void newFace();
    Face &getLastFace();
    void drawBorder();
    void drawColor();
    void translate(double dx, double dy, double dz);
    void rotate(double ux, double uy, double uz, double angle);
    void scale(double fx, double fy, double fz);
    void addDependency(Mesh &m);
    void writeToDisk(string filename);

    void setAmbientColor(double red, double green, double blue, double alpha);
    void setSpecularColor(double red, double green, double blue, double alpha);
    void setDiffuseColor(double red, double green, double blue, double alpha);
    void setShininess(double value);
};

// Methods implementation

Vector Vector::operator+(Vector &vect)
{
    return Vector(x + vect.x, y + vect.y, z + vect.z);
}

Vector Vector::operator-(Vector &vect)
{
    return Vector(x - vect.x, y - vect.y, z - vect.z);
}

Vector Vector::operator-() const
{
    return Vector(-x, -y, -z);
}

// Scala product
Vector Vector::operator*(double factor)
{
    return Vector(factor * x, factor * y, factor * z);
}
Vector Vector::operator*(double factor) const
{
    return Vector(factor * x, factor * y, factor * z);
}

// Scala division
Vector Vector::operator/(double divisor)
{
    return Vector(x / divisor, y / divisor, z / divisor);
}

// Dot product
double Vector::operator*(Vector &vect)
{
    return x * vect.x + y * vect.y + z * vect.z;
}

// Cross product
Vector Vector::operator|(Vector &vect)
{
    double a1 = x, a2 = y, a3 = z, b1 = vect.x, b2 = vect.y, b3 = vect.z;
    return Vector(b3 * a2 - b2 * a3, b1 * a3 - b3 * a1, b2 * a1 - b1 * a2);
}

// Length
double Vector::len()
{
    return sqrt(x * x + y * y + z * z);
}

// Normalize
Vector Vector::normalize()
{
    if (len() < ZERO)
        return *this;
    double l = len();
    return Vector(x / l, y / l, z / l);
}

// Print
void Vector::print()
{
    cout << setprecision(2) << fixed << x << "\t" << y << "\t" << z << endl;
}

bool Vector::isZero()
{
    return len() < ZERO;
};

void Face::setMainVertices(vector<Point *> &mainVertices)
{
    mainVertice = &mainVertices;
}

// void Face::setColor(double red, double green, double blue, double alpha)
// {
//     cr = red;
//     cb = blue;
//     cg = green;
//     ca = alpha;
// }

void Face::setAmbientColor(double red, double green, double blue, double alpha)
{
    ca[0] = red;
    ca[1] = green;
    ca[2] = blue;
    ca[3] = alpha;
}

void Face::setDiffuseColor(double red, double green, double blue, double alpha)
{
    cd[0] = red;
    cd[1] = green;
    cd[2] = blue;
    cd[3] = alpha;
}

void Face::setSpecularColor(double red, double green, double blue, double alpha)
{
    cs[0] = red;
    cs[1] = green;
    cs[2] = blue;
    cs[3] = alpha;
}

void Face::setShininess(double value)
{
    shininess = value;
}

void Face::addPoint(int index)
{
    vertice.push_back(index);
}

void Face::recalNorm()
{
    double _x = 0, _y = 0, _z = 0;
    for (int i = 0; i < (int)(vertice.size()); i++)
    {
        Point &p = *(mainVertice->at(vertice[i]));
        Point &pnext = *(mainVertice->at((vertice[(i + 1) % vertice.size()])));
        _x += (p.y - pnext.y) * (p.z + pnext.z);
        _y += (p.z - pnext.z) * (p.x + pnext.x);
        _z += (p.x - pnext.x) * (p.y + pnext.y);
    }
    norm = Vector(_x, _y, _z).normalize();
    if (mesh->reverse)
        norm = -norm;
}

void Face::drawBorder()
{
    glColor4d(ca[0], ca[1], ca[2], ca[3]);
    glBegin(GL_LINE_STRIP);
    {
        // for (vector<int>::iterator i = vertice.begin(); i != vertice.end(); i++)
        for (int i = 0; i < (int)(vertice.size()); i++)
        {
            Point &p = *(mainVertice->at(vertice[i]));
            // p.print();
            glVertex3d(p.x, p.y, p.z);
        }
        Point &p = *(mainVertice->at(vertice[0]));
        // p.print();
        // cout << endl;
        glVertex3d(p.x, p.y, p.z);
    }
    glEnd();
}

void Face::drawColor()
{
    recalNorm();
    glBegin(GL_POLYGON);
    {
        glNormal3d(norm.x, norm.y, norm.z);
        glColor4fv(ca);
        glMaterialfv(GL_FRONT, GL_AMBIENT, ca);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, cd);
        glMaterialfv(GL_FRONT, GL_SPECULAR, cs);
        glMaterialf(GL_FRONT, GL_SHININESS, shininess);
        // for (vector<int>::iterator i = vertice.begin(); i != vertice.end(); i++)
        for (int i = 0; i < (int)(vertice.size()); i++)
        {
            Point &p = *(mainVertice->at(vertice[i]));
            glVertex3d(p.x, p.y, p.z);
        }
        Point &p = *(mainVertice->at(vertice[0]));
        glVertex3d(p.x, p.y, p.z);
    }
    glEnd();
}

void Face::setMesh(Mesh *m)
{
    mesh = m;
}

void Mesh::addPoint(Point &point)
{
    mainVertice.push_back(new Point(point));
}

void Mesh::addFace(Face &_face)
{
    face.push_back(&_face);
    _face.setMesh(this);
}

void Mesh::newFace()
{
    face.push_back(new Face(&mainVertice));
    face[face.size() - 1]->setMesh(this);
}

Face &Mesh::getLastFace()
{
    return *(face[face.size() - 1]);
}

void Mesh::drawBorder()
{
    // for (vector<Face *>::iterator f = face.begin(); f != face.end(); f++)
    for (int i = 0; i < (int)(face.size()); i++)
        face[i]->drawBorder();
}

void Mesh::drawColor()
{
    for (int i = 0; i < (int)(face.size()); i++)
        face[i]->drawColor();
}

void Mesh::addDependency(Mesh &m)
{
    m.dpO = this;
    dp.push_back(&m);
}

void Mesh::translate(double dx, double dy, double dz)
{
    Vector tranV(dx, dy, dz);
    for (int i = 0; i < (int)(mainVertice.size()); i++)
    {
        *mainVertice[i] = *(mainVertice[i]) + tranV;
    }

    // recal normal vector
    for (int i = 0; i < (int)(face.size()); i++)
        face[i]->recalNorm();

    for (int i = 0; i < (int)(dp.size()); i++)
        dp[i]->translate(dx, dy, dz);
}

void Mesh::rotate(double ux, double uy, double uz, double angle)
{
    double c = cos(angle);
    double s = sin(angle);

    double r11 = c + (1 - c) * ux * ux;
    double r12 = (1 - c) * uy * ux - s * uz;
    double r13 = (1 - c) * uz * ux + s * uy;

    double r21 = (1 - c) * ux * uy + s * uz;
    double r22 = c + (1 - c) * uy * uy;
    double r23 = (1 - c) * uz * uy - s * ux;

    double r31 = (1 - c) * ux * uz - s * uy;
    double r32 = (1 - c) * uy * uz + s * ux;
    double r33 = c + (1 - c) * uz * uz;

    for (int i = 0; i < (int)(mainVertice.size()); i++)
    {
        Vector &v = *mainVertice[i];
        v = Point(
            v.x * r11 + v.y * r12 + v.z * r13,
            v.x * r21 + v.y * r22 + v.z * r23,
            v.x * r31 + v.y * r32 + v.z * r33);
    }

    // recal normal vector
    for (int i = 0; i < (int)face.size(); i++)
        face[i]->recalNorm();

    for (int i = 0; i < (int)dp.size(); i++)
        dp[i]->rotate(ux, uy, uz, angle);
}

void Mesh::scale(double fx, double fy, double fz)
{
    if (fx * fy * fz < 0)
        reverse = !reverse;

    for (int i = 0; i < (int)mainVertice.size(); i++)
    {
        Vector &v = *mainVertice[i];
        v = Point(fx * v.x, fy * v.y, fz * v.z);
    }

    // recal normal vector
    for (int i = 0; i < (int)(face.size()); i++)
        face[i]->recalNorm();
    for (int i = 0; i < (int)(dp.size()); i++)
        dp[i]->scale(fx, fy, fz);
}

void Mesh::setAmbientColor(double red, double green, double blue, double alpha)
{

    // for (vector<Face *>::iterator f = face.begin(); f != face.end(); f++)
    for (int i = 0; i < (int)(face.size()); i++)
        face[i]->setAmbientColor(red, green, blue, alpha);
};
void Mesh::setSpecularColor(double red, double green, double blue, double alpha)
{

    // for (vector<Face *>::iterator f = face.begin(); f != face.end(); f++)
    for (int i = 0; i < (int)(face.size()); i++)
        face[i]->setSpecularColor(red, green, blue, alpha);
};
;
void Mesh::setDiffuseColor(double red, double green, double blue, double alpha)
{

    // for (vector<Face *>::iterator f = face.begin(); f != face.end(); f++)
    for (int i = 0; i < (int)(face.size()); i++)
        face[i]->setDiffuseColor(red, green, blue, alpha);
};
;
void Mesh::setShininess(double value)
{

    // for (vector<Face *>::iterator f = face.begin(); f != face.end(); f++)
    for (int i = 0; i < (int)(face.size()); i++)
        face[i]->setShininess(value);
};

Mesh &createHolder(double bottomX, double bottomY, double bottomZ,
                   double lx, double lz, double height, double inR,
                   int acc)
{
    Point *fdr = new Point(lx / 2, 0, -lz / 2);
    Point *fdl = new Point(-lx / 2, 0, -lz / 2);
    Point *bdl = new Point(-lx / 2, 0, lz / 2);
    Point *bdr = new Point(lx / 2, 0, lz / 2);

    Point *fur = new Point(lx / 2, height - lx / 2, -lz / 2);
    Point *ful = new Point(-lx / 2, height - lx / 2, -lz / 2);
    Point *bul = new Point(-lx / 2, height - lx / 2, lz / 2);
    Point *bur = new Point(lx / 2, height - lx / 2, lz / 2);

    Point *fir = new Point(inR, height - lx / 2, -lz / 2);
    Point *fil = new Point(-inR, height - lx / 2, -lz / 2);
    Point *bil = new Point(-inR, height - lx / 2, lz / 2);
    Point *bir = new Point(inR, height - lx / 2, lz / 2);

    Point *fjr = new Point(inR, 0, -lz / 2);
    Point *fjl = new Point(-inR, 0, -lz / 2);
    Point *bjl = new Point(-inR, 0, lz / 2);
    Point *bjr = new Point(inR, 0, lz / 2);

    int &N = acc;
    map<void *, int> pMap;
    Mesh &mesh = *(new Mesh());
    int index = 0;

    {
        mesh.addPoint(*fdr);
        pMap[fdr] = index++;

        mesh.addPoint(*fdl);
        pMap[fdl] = index++;

        mesh.addPoint(*bdl);
        pMap[bdl] = index++;

        mesh.addPoint(*bdr);
        pMap[bdr] = index++;

        mesh.addPoint(*ful);
        pMap[ful] = index++;

        mesh.addPoint(*fur);
        pMap[fur] = index++;

        mesh.addPoint(*bul);
        pMap[bul] = index++;

        mesh.addPoint(*bur);
        pMap[bur] = index++;

        mesh.addPoint(*fil);
        pMap[fil] = index++;

        mesh.addPoint(*fir);
        pMap[fir] = index++;

        mesh.addPoint(*bil);
        pMap[bil] = index++;

        mesh.addPoint(*bir);
        pMap[bir] = index++;

        mesh.addPoint(*fjl);
        pMap[fjl] = index++;

        mesh.addPoint(*fjr);
        pMap[fjr] = index++;

        mesh.addPoint(*bjl);
        pMap[bjl] = index++;

        mesh.addPoint(*bjr);
        pMap[bjr] = index++;
    }

    vector<Point *> fi, fo, bi, bo, fd, bd;
    for (int i = 0; i < (int)N; i++)
    {
        double angle = 2.0 * M_PI / N * i;
        Point *pf = new Point(inR * cos(angle), height - lx / 2 + inR * sin(angle), -lz / 2);
        Point *pb = new Point(inR * cos(angle), height - lx / 2 + inR * sin(angle), lz / 2);
        fi.push_back(pf);
        bi.push_back(pb);
        mesh.addPoint(*pf);
        mesh.addPoint(*pb);
        pMap[pf] = index++;
        pMap[pb] = index++;
    }
    fi.push_back(fi[0]);
    bi.push_back(bi[0]);

    for (int i = N / 2; i < (int)N + 1; i++)
    {
        double angle = 2.0 * M_PI / N * i;
        Point *pf = new Point(inR * cos(angle), 0, -lz / 2);
        Point *pb = new Point(inR * cos(angle), 0, lz / 2);
        fd.push_back(pf);
        bd.push_back(pb);
        mesh.addPoint(*pf);
        mesh.addPoint(*pb);
        pMap[pf] = index++;
        pMap[pb] = index++;
    }

    double outR = lx / 2;
    for (int i = 0; i < (int)N / 2 + 1; i++)
    {
        double angle = 2.0 * M_PI / N * i;
        Point *pf = new Point(outR * cos(angle), height - lx / 2 + outR * sin(angle), -lz / 2);
        Point *pb = new Point(outR * cos(angle), height - lx / 2 + outR * sin(angle), lz / 2);
        fo.push_back(pf);
        bo.push_back(pb);
        mesh.addPoint(*pf);
        mesh.addPoint(*pb);
        pMap[pf] = index++;
        pMap[pb] = index++;
    }
    fo.push_back(fo[0]);
    bo.push_back(bo[0]);

    // Bottom face
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[fdr]);
        face.addPoint(pMap[bdr]);
        face.addPoint(pMap[bdl]);
        face.addPoint(pMap[fdl]);
    }

    // Left face
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[fdl]);
        face.addPoint(pMap[bdl]);
        face.addPoint(pMap[bul]);
        face.addPoint(pMap[ful]);
    }

    // Right face
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[fdr]);
        face.addPoint(pMap[fur]);
        face.addPoint(pMap[bur]);
        face.addPoint(pMap[bdr]);
    }

    // Front face 1
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[fdl]);
        face.addPoint(pMap[ful]);
        face.addPoint(pMap[fil]);
        face.addPoint(pMap[fjl]);
    }

    // Front face 2
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[fdr]);
        face.addPoint(pMap[fjr]);
        face.addPoint(pMap[fir]);
        face.addPoint(pMap[fur]);
    }

    // Back face 1
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[bdl]);
        face.addPoint(pMap[bjl]);
        face.addPoint(pMap[bil]);
        face.addPoint(pMap[bul]);
    }

    // Back face 2
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[bdr]);
        face.addPoint(pMap[bur]);
        face.addPoint(pMap[bir]);
        face.addPoint(pMap[bjr]);
    }
    // Outside face
    for (int i = 0; i < (int)N / 2; i++)
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[fo[i]]);
        face.addPoint(pMap[fo[i + 1]]);
        face.addPoint(pMap[bo[i + 1]]);
        face.addPoint(pMap[bo[i]]);
    }

    // Inside face
    for (int i = 0; i < (int)N; i++)
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[fi[i]]);
        face.addPoint(pMap[bi[i]]);
        face.addPoint(pMap[bi[(i + 1) % N]]);
        face.addPoint(pMap[fi[(i + 1) % N]]);
    }

    // Front face
    for (int i = 0; i < (int)N / 2; i++)
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[fi[i]]);
        face.addPoint(pMap[fi[i + 1]]);
        face.addPoint(pMap[fo[i + 1]]);
        face.addPoint(pMap[fo[i]]);
    }

    for (int i = N / 2; i < (int)N; i++)
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[fi[i]]);
        face.addPoint(pMap[fi[i + 1]]);
        face.addPoint(pMap[fd[i - N / 2 + 1]]);
        face.addPoint(pMap[fd[i - N / 2]]);
    }

    // Back face
    for (int i = 0; i < (int)N / 2; i++)
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[bi[i]]);
        face.addPoint(pMap[bo[i]]);
        face.addPoint(pMap[bo[i + 1]]);
        face.addPoint(pMap[bi[i + 1]]);
    }

    for (int i = N / 2; i < (int)N - 1; i++)
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[bi[i]]);
        face.addPoint(pMap[bd[i - N / 2]]);
        face.addPoint(pMap[bd[i - N / 2 + 1]]);
        face.addPoint(pMap[bi[i + 1]]);
    }
    mesh.translate(bottomX, bottomY, bottomZ);

    // Clean up
    delete fdr;
    delete fdl;
    delete bdl;
    delete bdr;

    delete fur;
    delete ful;
    delete bul;
    delete bur;

    delete fir;
    delete fil;
    delete bil;
    delete bir;

    delete fjr;
    delete fjl;
    delete bjl;
    delete bjr;

    pMap.clear();
    for (int i = 0; i < (int)fi.size() - 1; i++)
        delete fi[i];
    for (int i = 0; i < (int)fo.size() - 1; i++)
        delete fo[i];
    for (int i = 0; i < (int)fd.size() - 1; i++)
        delete fd[i];
    for (int i = 0; i < (int)bi.size() - 1; i++)
        delete bi[i];
    for (int i = 0; i < (int)bo.size() - 1; i++)
        delete bo[i];
    for (int i = 0; i < (int)bd.size() - 1; i++)
        delete bd[i];
    fi.clear();
    fo.clear();
    fd.clear();
    bi.clear();
    bo.clear();
    bd.clear();

    return mesh;
}

Mesh &createCylinder(double bottomX, double bottomY, double bottomZ,
                     double radius, double height,
                     int acc)
{

    int &N = acc;
    Mesh &mesh = *(new Mesh("empty.dat"));
    Point *topCenter = new Point(0, height, 0);
    Point *bottomCenter = new Point(0, 0, 0);
    mesh.addPoint(*topCenter);
    mesh.addPoint(*bottomCenter);

    map<void *, int> pMap;
    pMap[topCenter] = 0;
    pMap[bottomCenter] = 1;
    int index = 2;
    Point *topPoint = new Point[N]();
    Point *bottomPoint = new Point[N]();

    for (int i = 0; i < N; i++)
    {
        double angle = 2.0 * M_PI / N * i;
        topPoint[i] = Point(radius * cos(angle), height, radius * sin(angle));
        bottomPoint[i] = Point(radius * cos(angle), 0, radius * sin(angle));
        mesh.addPoint(topPoint[i]);
        mesh.addPoint(bottomPoint[i]);
        pMap[topPoint + i] = index++;
        pMap[bottomPoint + i] = index++;
    }

    for (int i = N - 1; i >= 0; i--)
    {
        // Top face
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[topPoint + (i + 1) % N]);
        face.addPoint(pMap[topPoint + i]);

        face.addPoint(0);
    }

    for (int i = 0; i < N; i++)
    {

        {
            // Bottom face
            mesh.newFace();
            Face &face = mesh.getLastFace();
            face.addPoint(pMap[bottomPoint + i]);
            face.addPoint(pMap[bottomPoint + (i + 1) % N]);
            face.addPoint(1);
        }
        {
            // Side face
            mesh.newFace();
            Face &face = mesh.getLastFace();
            face.addPoint(pMap[topPoint + i]);
            face.addPoint(pMap[topPoint + (i + 1) % N]);
            face.addPoint(pMap[bottomPoint + (i + 1) % N]);
            face.addPoint(pMap[bottomPoint + i]);
        }
    }

    mesh.translate(bottomX, bottomY, bottomZ);

    // Clean up
    delete[] topPoint;
    delete[] bottomPoint;
    delete topCenter;
    delete bottomCenter;

    return mesh;
}

Mesh &createRectangular(double bottomX, double bottomY, double bottomZ,
                        double rx, double ry, double rz)
{
    Point *bsw = new Point(-rx / 2, 0, -rz / 2);
    Point *bse = new Point(rx / 2, 0, -rz / 2);
    Point *bnw = new Point(-rx / 2, 0, rz / 2);
    Point *bne = new Point(rx / 2, 0, rz / 2);

    Point *usw = new Point(-rx / 2, ry, -rz / 2);
    Point *use = new Point(rx / 2, ry, -rz / 2);
    Point *unw = new Point(-rx / 2, ry, rz / 2);
    Point *une = new Point(rx / 2, ry, rz / 2);

    map<void *, int> pMap;
    Mesh &mesh = *(new Mesh());
    int index = 0;

    mesh.addPoint(*bsw);
    pMap[bsw] = index++;

    mesh.addPoint(*bse);
    pMap[bse] = index++;

    mesh.addPoint(*bnw);
    pMap[bnw] = index++;

    mesh.addPoint(*bne);
    pMap[bne] = index++;

    mesh.addPoint(*usw);
    pMap[usw] = index++;

    mesh.addPoint(*use);
    pMap[use] = index++;

    mesh.addPoint(*unw);
    pMap[unw] = index++;

    mesh.addPoint(*une);
    pMap[une] = index++;
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[usw]);
        face.addPoint(pMap[use]);
        face.addPoint(pMap[bse]);
        face.addPoint(pMap[bsw]);
    }
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[use]);
        face.addPoint(pMap[une]);
        face.addPoint(pMap[bne]);
        face.addPoint(pMap[bse]);
    }
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[une]);
        face.addPoint(pMap[unw]);
        face.addPoint(pMap[bnw]);
        face.addPoint(pMap[bne]);
    }
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[unw]);
        face.addPoint(pMap[usw]);
        face.addPoint(pMap[bsw]);
        face.addPoint(pMap[bnw]);
    }
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[unw]);
        face.addPoint(pMap[une]);
        face.addPoint(pMap[use]);
        face.addPoint(pMap[usw]);
    }
    {
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[bnw]);
        face.addPoint(pMap[bsw]);
        face.addPoint(pMap[bse]);
        face.addPoint(pMap[bne]);
    }
    mesh.translate(bottomX, bottomY, bottomZ);
    return mesh;
}

Mesh &createDuoCylinder(double bottomX, double bottomY, double bottomZ,
                        double h1, double hs, double h2, double rd, double ru,
                        int acc)
{
    int &N = acc;
    Mesh &mesh = *(new Mesh("empty.dat"));
    Point *bottomCenter = new Point(0, 0, 0);
    Point *topCenter = new Point(0, h1 + hs + h2, 0);
    mesh.addPoint(*bottomCenter);
    mesh.addPoint(*topCenter);
    map<void *, int> pMap;
    pMap[bottomCenter] = 0;
    pMap[topCenter] = 1;
    int index = 2;
    Point *l0 = new Point[N]();
    Point *l1 = new Point[N]();
    Point *l2 = new Point[N]();
    Point *l3 = new Point[N]();
    for (int i = 0; i < N; i++)
    {
        double angle = 2.0 * M_PI / N * i;
        l0[i] = Point(rd * cos(angle), 0, rd * sin(angle));
        l1[i] = Point(rd * cos(angle), h1, rd * sin(angle));
        l2[i] = Point(ru * cos(angle), h1 + hs, ru * sin(angle));
        l3[i] = Point(ru * cos(angle), h1 + hs + h2, ru * sin(angle));
        mesh.addPoint(l0[i]);
        mesh.addPoint(l1[i]);
        mesh.addPoint(l2[i]);
        mesh.addPoint(l3[i]);
        pMap[l0 + i] = index++;
        pMap[l1 + i] = index++;
        pMap[l2 + i] = index++;
        pMap[l3 + i] = index++;
    }

    for (int i = 0; i < N; i++)
    {
        {
            // Top face
            mesh.newFace();
            Face &face = mesh.getLastFace();
            face.addPoint(pMap[l3 + i]);
            face.addPoint(pMap[topCenter]);
            face.addPoint(pMap[l3 + (i + 1) % N]);
        }
        {
            // Bottom face
            mesh.newFace();
            Face &face = mesh.getLastFace();
            face.addPoint(pMap[l0 + (i + 1) % N]);
            face.addPoint(pMap[bottomCenter]);
            face.addPoint(pMap[l0 + i]);
        }
        {
            // Bottom side face
            mesh.newFace();
            Face &face = mesh.getLastFace();
            face.addPoint(pMap[l1 + i]);
            face.addPoint(pMap[l1 + (i + 1) % N]);
            face.addPoint(pMap[l0 + (i + 1) % N]);
            face.addPoint(pMap[l0 + i]);
        }
        {
            // Middle side face
            mesh.newFace();
            Face &face = mesh.getLastFace();
            face.addPoint(pMap[l2 + i]);
            face.addPoint(pMap[l2 + (i + 1) % N]);
            face.addPoint(pMap[l1 + (i + 1) % N]);
            face.addPoint(pMap[l1 + i]);
        }
        {
            // Top side face
            mesh.newFace();
            Face &face = mesh.getLastFace();
            face.addPoint(pMap[l3 + i]);
            face.addPoint(pMap[l3 + (i + 1) % N]);
            face.addPoint(pMap[l2 + (i + 1) % N]);
            face.addPoint(pMap[l2 + i]);
        }
    }

    mesh.translate(bottomX, bottomY, bottomZ);
    // Clean up
    delete[] l0;
    delete[] l1;
    delete[] l2;
    delete[] l3;
    // delete topCenter;
    // delete bottomCenter;
    return mesh;
}

Mesh &createSpeacialOval(double bottomX, double bottomY, double bottomZ,
                         double ri, double ro, double len, double height,
                         int acc)
{
    int &N = acc;
    Mesh &mesh = *(new Mesh("empty.dat"));
    map<void *, int> pMap;
    Point *_ui = new Point[N + 2]();
    Point *_uo = new Point[N + 2]();
    Point *_di = new Point[N + 2]();
    Point *_do = new Point[N + 2]();
    int index = 0;
    for (int i = 0; i <= N / 2; i++)
    {
        double angle = 2.0 * M_PI / N * i;
        _ui[i] = Point(ri * cos(angle), height, ri * sin(angle) + len / 2);
        mesh.addPoint(_ui[i]);
        pMap[_ui + i] = index++;
        _uo[i] = Point(ro * cos(angle), height, ro * sin(angle) + len / 2);
        mesh.addPoint(_uo[i]);
        pMap[_uo + i] = index++;
        _di[i] = Point(ri * cos(angle), 0, ri * sin(angle) + len / 2);
        mesh.addPoint(_di[i]);
        pMap[_di + i] = index++;
        _do[i] = Point(ro * cos(angle), 0, ro * sin(angle) + len / 2);
        mesh.addPoint(_do[i]);
        pMap[_do + i] = index++;
    }
    for (int i = N / 2 + 1; i <= N + 1; i++)
    {
        double angle = 2.0 * M_PI / N * (i - 1);
        _ui[i] = Point(ri * cos(angle), height, ri * sin(angle) - len / 2);
        mesh.addPoint(_ui[i]);
        pMap[_ui + i] = index++;
        _uo[i] = Point(ro * cos(angle), height, ro * sin(angle) - len / 2);
        mesh.addPoint(_uo[i]);
        pMap[_uo + i] = index++;
        _di[i] = Point(ri * cos(angle), 0, ri * sin(angle) - len / 2);
        mesh.addPoint(_di[i]);
        pMap[_di + i] = index++;
        _do[i] = Point(ro * cos(angle), 0, ro * sin(angle) - len / 2);
        mesh.addPoint(_do[i]);
        pMap[_do + i] = index++;
    }

    for (int i = 0; i < N + 2; i++)
    {
        // Top face
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[_ui + i]);
        face.addPoint(pMap[_ui + (i + 1) % (N + 2)]);
        face.addPoint(pMap[_uo + (i + 1) % (N + 2)]);
        face.addPoint(pMap[_uo + i]);
    }
    for (int i = 0; i < N + 2; i++)
    {
        // Bottom face
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[_di + i]);
        face.addPoint(pMap[_do + i]);
        face.addPoint(pMap[_do + (i + 1) % (N + 2)]);
        face.addPoint(pMap[_di + (i + 1) % (N + 2)]);
    }
    for (int i = 0; i < N + 2; i++)
    {
        // Inside face
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[_ui + i]);
        face.addPoint(pMap[_di + i]);
        face.addPoint(pMap[_di + (i + 1) % (N + 2)]);
        face.addPoint(pMap[_ui + (i + 1) % (N + 2)]);
    }
    for (int i = 0; i < N + 2; i++)
    {
        // Outside face
        mesh.newFace();
        Face &face = mesh.getLastFace();
        face.addPoint(pMap[_uo + i]);
        face.addPoint(pMap[_uo + (i + 1) % (N + 2)]);
        face.addPoint(pMap[_do + (i + 1) % (N + 2)]);
        face.addPoint(pMap[_do + i]);
    }

    mesh.translate(bottomX, bottomY, bottomZ);
    delete[] _uo;
    delete[] _ui;
    delete[] _do;
    delete[] _di;
    return mesh;
}

//=-----------

#define N 120

#define AMB_COPPER 0.19125, 0.0735, 0.0225, 1
#define DIF_COPPER 0.7038, 0.37048, 0.0828, 1
#define SPE_COPPER 0.256777, 0.137622, 0.086014, 1
#define SHI_COPPER 0.1

#define AMB_GOLD 0.24725, 0.1995, 0.0745, 1
#define DIF_GOLD 0.75164, 0.60648, 0.22648, 1
#define SPE_GOLD 0.628281, 0.555802, 0.366065, 1
#define SHI_GOLD 0.4

#define AMB_SILVER 0.19225, 0.19225, 0.19225, 1
#define DIF_SILVER 0.50754, 0.50754, 0.50754, 1
#define SPE_SILVER 0.508273, 0.508273, 0.508273, 1
#define SHI_SILVER 0.4

#define AMB_OBSIDIAN 0.05375, 0.05, 0.06625, 1
#define DIF_OBSIDIAN 0.18275, 0.17, 0.22525, 1
#define SPE_OBSIDIAN 0.532741, 0.528634, 0.546435, 1
#define SHI_OBSIDIAN 0.6

#define AMB_CHROME 0.25, 0.25, 0.25, 1
#define DIF_CHROME 0.4, 0.4, 0.4, 1
#define SPE_CHROME 0.774597, 0.774597, 0.774597, 1
#define SHI_CHROME 0.6

#define AMB_CYAN_PLASTIC 0.0, 0.1, 0.06, 1
#define DIF_CYAN_PLASTIC 0.0, 0.50980392, 0.50980392, 1
#define SPE_CYAN_PLASTIC 0.50196078, 0.50196078, 0.50196078, 1
#define SHI_CYAN_PLASTIC 0.25

#define AMB_RUBY 0.1745, 0.01175, 0.01175, 1
#define DIF_RUBY 0.61424, 0.04136, 0.04136, 1
#define SPE_RUBY 0.727811, 0.626959, 0.626959, 1
#define SHI_RUBY 0.6

#define AMB_TURQUOISE 0.1, 0.18725, 0.1745, 1
#define DIF_TURQUOISE 0.396, 0.74151, 0.69102, 1
#define SPE_TURQUOISE 0.297254, 0.30829, 0.306678, 1
#define SHI_TURQUOISE 0.1

#define AMB_GREEN_PLASTIC 0.0, 0.0, 0.0, 1
#define DIF_GREEN_PLASTIC 0.1, 0.75, 0.1, 1
#define SPE_GREEN_PLASTIC 0.45, 0.55, 0.45, 1
#define SHI_GREEN_PLASTIC .25

#define AMB_BLUE_PLASTIC 0.0, 0.0, 0.0, 1
#define DIF_BLUE_PLASTIC 0.1, 0.1, 0.65, 1
#define SPE_BLUE_PLASTIC 0.45, 0.55, 0.45, 1
#define SHI_BLUE_PLASTIC .25

bool bWireFrame = false;
GLfloat viewX = 30 * cos(M_PI / 4), viewY = 10, viewZ = 30 * sin(M_PI / 4);
GLfloat pistonHeight = 0;
GLfloat D = 30;
GLfloat angle = M_PI / 4;
GLfloat angle1 = 0;
GLfloat angle2 = 0;
GLfloat angle3 = 0;
double vx = 0, vy = 0, vz = 0;
double lx = 1, lz = 0;
bool light1 = false;

double dx = 0;

bool autoMove = false;

double ran()
{
    return (double(rand()) / RAND_MAX) * 2 - 1;
}

vector<Mesh *> mesh;

// Draw the floor
void drawFloor()
{
    glDisable(GL_LIGHTING);
    // // Begin drawing the floor
    // glBegin(GL_QUADS);

    // glColor4f(1, 1, 1, alpha);
    // glNormal3f(0, 1, 0);
    // glVertex3f(-width / 2, 0, length / 2);
    // glVertex3f(-width / 2, 0, -length / 2);
    // glVertex3f(width / 2, 0, -length / 2);
    // glVertex3f(width / 2, 0, length / 2);

    // glEnd();

    // Coordinates
    glBegin(GL_LINES);
    {
        glColor3f(1, 0, 0);
        glVertex3f(100, 0, 0);
        glVertex3f(0, 0, 0);
        glColor3f(0, 1, 0);
        glVertex3f(0, 100, 0);
        glVertex3f(0, 0, 0);

        glColor3f(0, 0, 1);
        glVertex3f(0, 0, 100);
        glVertex3f(0, 0, 0);
    }
    glEnd();

    double angle = M_PI / 6;
    double alpha = 0.5;
    double r = 2;
    for (double x = -50; x < 50; x += r * sqrt(3))
        for (double y = -50; y < 50; y += 3 * r)
            for (int i = 0; i < 3; i++)
            {
                glBegin(GL_POLYGON);
                {
                    glColor4f(i < 1 ? 1 : 0, i < 2 ? 1 : 0, 1, alpha);
                    glVertex3f(x, 0, y);
                    angle -= M_PI / 3;
                    for (int i = 0; i < 3; i++)
                    {
                        angle += M_PI / 3;
                        glVertex3f(x + cos(angle) * r, 0, y + sin(angle) * r);
                    }
                }
                glEnd();
            }

    for (double x = -50 - r * sqrt(3) / 2; x < 50; x += r * sqrt(3))
        for (double y = -50 - 3 * r / 2; y < 50; y += 3 * r)
            for (int i = 0; i < 3; i++)
            {
                glBegin(GL_POLYGON);
                {
                    glColor4f(i < 1 ? 1 : 0, i < 2 ? 1 : 0, 1, alpha);
                    glVertex3f(x, 0, y);
                    angle -= M_PI / 3;
                    for (int i = 0; i < 3; i++)
                    {
                        angle += M_PI / 3;
                        glVertex3f(x + cos(angle) * r, 0, y + sin(angle) * r);
                    }
                }
                glEnd();
            }

    glEnable(GL_LIGHTING);
}

void drawColor()
{
    glPushMatrix();
    for (int i = 0; i < (int)mesh.size(); i++)
        mesh[i]->drawColor();
    glPopMatrix();
}

void drawWireFrame()
{
    glDisable(GL_LIGHTING);
    glPushMatrix();
    for (int i = 0; i < (int)mesh.size(); i++)
        mesh[i]->drawBorder();
    glPopMatrix();
    glEnable(GL_LIGHTING);
}

void display()
{
    // Clear the stencil buffer
    glClearStencil(0);
    // Clear depth
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(viewX, viewY, viewZ, 0, 0, 0, 0, 1, 0);
    // Set light position
    GLfloat lightPosition[] = {15, 20, 15, 1};
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    glEnable(GL_LIGHT0);
    if (light1)
    {
        GLfloat lightPosition1[] = {-15, 3, 15, 1};
        glLightfv(GL_LIGHT1, GL_POSITION, lightPosition1);
        glEnable(GL_LIGHT1);
    }
    else
    {
        glDisable(GL_LIGHT1);
    }

    // Draw
    glPushMatrix();
    if (!bWireFrame)
        drawColor();
    else
        drawWireFrame();
    glPopMatrix();

    // Enable the stencil buffer
    glEnable(GL_STENCIL_TEST);
    // Disable drawing colors
    glColorMask(0, 0, 0, 0);
    // Disable depth testing
    glDisable(GL_DEPTH_TEST);
    // Make the stencil test always pass
    glStencilFunc(GL_ALWAYS, 1, 1);
    // Make pixels in the stencil buffer be set to 1 when the stencil test passes
    glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
    // Set all of the pixels covered by the floor to be 1 in the stencil buffer
    drawFloor();

    // Enable drawing colors to the screen
    glColorMask(1, 1, 1, 1);
    // Enable depth testing
    glEnable(GL_DEPTH_TEST);
    // Make the stencil test pass only when the pixel is 1 in the stencil buffer
    glStencilFunc(GL_EQUAL, 1, 1);
    // Make the stencil buffer not change
    glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);

    //Draw the cube, reflected vertically, at all pixels where the stencil
    //buffer is 1
    glPushMatrix();
    glScalef(1, -1, 1);
    // glTranslatef(0, 5.0f, 0);
    if (!bWireFrame)
        drawColor();
    else
        drawWireFrame();
    glPopMatrix();

    // Disable using the stencil buffer
    glDisable(GL_STENCIL_TEST);

    // Blend the floor onto the screen
    glEnable(GL_BLEND);
    drawFloor();
    glDisable(GL_BLEND);

    glutSwapBuffers();
}

void reshape(int width, int height)
{
    int size = width < height ? width : height;
    glViewport(0, 0, size, size);
}

void myKeyboard(unsigned char key, int x, int y)
{
    static GLfloat angle1Step = 0.03;
    static GLfloat angle2Step = 0.03;
    static GLfloat angle3Step = 0.01;
    static GLfloat pistonStep = 0.1;
    static GLfloat pistonMaxDiff = 4;

    Vector lookVector(viewX, viewY, viewZ);
    double l = lookVector.len();
    lookVector = lookVector.normalize();
    lookVector = lookVector * 0.5;
    int a = x;
    a = y;
    a += a;

    switch (key)
    {
    case '1':
        angle1 += angle1Step;
        mesh[0]->rotate(0, 1, 0, angle1Step);

        if (angle1 > 2 * M_PI)
            angle1 -= 2 * M_PI;
        break;
    case '2':
        angle1 -= angle1Step;
        mesh[0]->rotate(0, 1, 0, -angle1Step);
        if (angle1 < 0)
            angle1 += 2 * M_PI;
        break;
    case '3':
        angle2 += angle2Step;
        mesh[2]->rotate(0, 1, 0, angle2Step);
        if (angle2 > 2 * M_PI)
            angle2 -= 2 * M_PI;
        break;
    case '4':
        angle2 -= angle2Step;
        mesh[2]->rotate(0, 1, 0, -angle2Step);
        if (angle2 < 0)
            angle2 += 2 * M_PI;
        break;
    case '5':
        pistonHeight += pistonStep;
        mesh[2]->translate(0, pistonStep, 0);
        while (pistonHeight > pistonMaxDiff)
        {
            pistonHeight -= pistonStep;
            mesh[2]->translate(0, -pistonStep, 0);
        }
        break;
    case '6':
        pistonHeight -= pistonStep;
        mesh[2]->translate(0, -pistonStep, 0);
        while (pistonHeight < 0)
        {
            pistonHeight += pistonStep;
            mesh[2]->translate(0, pistonStep, 0);
        }
        break;
    case '7':
        for (int i = 0; i < 7; i++)
        {
            angle3 += angle3Step;
            mesh[4]->rotate(0, 1, 0, angle3Step);
            lx = cos(angle1 + angle2);
            lz = -sin(angle1 + angle2);
            double diff = 4 * sin(angle3) - dx;
            dx += diff;
            mesh[11]->translate(diff * lx, 0, diff * lz);
            if (angle3 > 2 * M_PI)
                angle3 -= 2 * M_PI;
        }
        break;
    case '8':
        for (int i = 0; i++ < 7;)
        {
            angle3 -= angle3Step;
            mesh[4]->rotate(0, 1, 0, -angle3Step);
            lx = cos(angle1 + angle2);
            lz = -sin(angle1 + angle2);
            double diff = 4 * sin(angle3) - dx;
            dx += diff;
            mesh[11]->translate(diff * lx, 0, diff * lz);
            if (angle3 < 0)
                angle3 += 2 * M_PI;
        }
        break;
    case 'w':
    case 'W':
        bWireFrame = !bWireFrame;
        break;
    case 'a':
    case 'A':
        autoMove = !autoMove;
        break;
    case 'd':
    case 'D':
        light1 = !light1;
        break;
    case '-':
        if (l < 100)
        {
            viewX += lookVector.x;
            viewY += lookVector.y;
            viewZ += lookVector.z;
        }

        break;
    case '+':
        if (l > 10)
        {
            viewX -= lookVector.x;
            viewY -= lookVector.y;
            viewZ -= lookVector.z;
        }

        break;
    }
    glutPostRedisplay();
}

void time(int)
{
    glutTimerFunc(50, time, 0);
    if (autoMove)
    {
        double angle3Step = 0.05;
        for (int i = 0; i < 2; i++)
        {
            angle3 += angle3Step;
            mesh[4]->rotate(0, 1, 0, angle3Step);
            lx = cos(angle1 + angle2);
            lz = -sin(angle1 + angle2);
            double diff = 4 * sin(angle3) - dx;
            dx += diff;
            mesh[11]->translate(diff * lx, 0, diff * lz);

            if (angle3 > 2 * M_PI)
                angle3 -= 2 * M_PI;
        }
    }
    glutPostRedisplay();
}

void mySpecialKeyboard(int key, int x, int y)
{
    int a = x;
    a = y;
    a += a;

    Vector lookVector(viewX, 0, viewZ);
    double dp = 0.05;
    double r = lookVector.len();
    switch (key)
    {
    case GLUT_KEY_UP:
        viewY += 0.5;
        while (viewY >= 60)
            viewY -= 0.5;
        break;
    case GLUT_KEY_DOWN:
        viewY -= 0.5;
        while (viewY <= 1)
            viewY += 0.5;
        break;
    case GLUT_KEY_RIGHT:
        angle += dp;
        viewX = r * cos(angle);
        viewZ = r * sin(angle);
        break;
    case GLUT_KEY_LEFT:
        angle -= dp;
        viewX = r * cos(angle);
        viewZ = r * sin(angle);
        break;
    default:
        break;
    }
    if (angle < 0)
        angle += 2 * M_PI;
    if (angle > 2 * M_PI)
        angle -= 2 * M_PI;
    glutPostRedisplay();
}

void init()
{
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_NORMALIZE);
    // glEnable(GL_COLOR_MATERIAL);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);

    // Setup light
    GLfloat lightAmbient[] = {0.4, .4, .4, 1};
    GLfloat lightDiffuse[] = {1, 1, 1, 1};
    GLfloat lightSpecular[] = {1, 1, 1, 1};

    GLfloat lightAmbient1[] = {0.4, 0.4, 0.1, 1};
    GLfloat lightDiffuse1[] = {0.3, 0.3, 0.3, 1};
    GLfloat lightSpecular1[] = {1, 1, 1, 1};

    glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);
    glLightfv(GL_LIGHT1, GL_AMBIENT, lightAmbient1);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, lightDiffuse1);
    glLightfv(GL_LIGHT1, GL_SPECULAR, lightSpecular1);

    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

    // glEnable(GL_COLOR_MATERIAL_FACE);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-2, 2, -2, 2, 1, 500); // Viewing volume

    // Timer
    glutTimerFunc(10, time, 0);

    glMatrixMode(GL_MODELVIEW);
    glClearColor(0, 0, 0, 1);

    // Draw object

    // Tripod 1
    mesh.push_back(&createCylinder(0, 0, 0, 4, 1, N));
    mesh[0]->setAmbientColor(AMB_RUBY);
    mesh[0]->setDiffuseColor(DIF_RUBY);
    mesh[0]->setSpecularColor(SPE_RUBY);
    mesh[0]->setShininess(100 * SHI_RUBY);

    // Tripod 2
    mesh.push_back(&createCylinder(0, 1, 0, 3, 5, N));
    mesh[1]->setAmbientColor(AMB_RUBY);
    mesh[1]->setDiffuseColor(DIF_RUBY);
    mesh[1]->setSpecularColor(SPE_RUBY);
    mesh[1]->setShininess(100 * SHI_RUBY);
    mesh[0]->addDependency(*mesh[1]);

    // Piston
    mesh.push_back(&createCylinder(0, 0, 0, 2.5, 6, N));
    mesh[2]->setAmbientColor(AMB_BLUE_PLASTIC);
    mesh[2]->setDiffuseColor(DIF_BLUE_PLASTIC);
    mesh[2]->setSpecularColor(SPE_BLUE_PLASTIC);
    mesh[2]->setShininess(100 * SHI_BLUE_PLASTIC);
    mesh[0]->addDependency(*mesh[2]);

    // Large panel
    mesh.push_back(&createRectangular(0, 0, 0, 40, 1, 15));
    mesh[3]->setAmbientColor(AMB_GREEN_PLASTIC);
    mesh[3]->setDiffuseColor(DIF_GREEN_PLASTIC);
    mesh[3]->setSpecularColor(SPE_GREEN_PLASTIC);
    mesh[3]->setShininess(100 * SHI_GREEN_PLASTIC);
    mesh[3]->translate(0, 6, 0);
    mesh[2]->addDependency(*mesh[3]);

    // Plate
    mesh.push_back(&createCylinder(0, 0, 0, 6, 1, N));
    mesh[4]->setAmbientColor(AMB_RUBY);
    mesh[4]->setDiffuseColor(DIF_RUBY);
    mesh[4]->setSpecularColor(SPE_RUBY);
    mesh[4]->setShininess(100 * SHI_RUBY);
    mesh[4]->translate(0, 7, 0);
    mesh[3]->addDependency(*mesh[4]);

    // Holder plates 1 2
    mesh.push_back(&createRectangular(0, 0, 0, 2, 0.1, 2));
    mesh.push_back(&createRectangular(0, 0, 0, 2, 0.1, 2));
    mesh[5]->setAmbientColor(AMB_OBSIDIAN);
    mesh[5]->setDiffuseColor(DIF_OBSIDIAN);
    mesh[5]->setSpecularColor(SPE_OBSIDIAN);
    mesh[5]->setShininess(100 * SHI_OBSIDIAN);
    mesh[5]->scale(1, 1, 2);
    mesh[5]->translate(10, 7, 0);
    mesh[3]->addDependency(*mesh[5]);

    mesh[6]->setAmbientColor(AMB_OBSIDIAN);
    mesh[6]->setDiffuseColor(DIF_OBSIDIAN);
    mesh[6]->setSpecularColor(SPE_OBSIDIAN);
    mesh[6]->setShininess(100 * SHI_OBSIDIAN);
    mesh[6]->scale(1, 1, 2);
    mesh[6]->translate(-10, 7, 0);
    mesh[3]->addDependency(*mesh[6]);

    // Holder 1 2
    mesh.push_back(&createHolder(0, 7.1, -10, 2, 2, 2.8, 0.4, N));
    mesh[7]->rotate(0, 1, 0, M_PI / 2);
    mesh[7]->setAmbientColor(AMB_OBSIDIAN);
    mesh[7]->setDiffuseColor(DIF_OBSIDIAN);
    mesh[7]->setSpecularColor(SPE_OBSIDIAN);
    mesh[7]->setShininess(100 * SHI_OBSIDIAN);
    mesh[5]->addDependency(*mesh[7]);

    mesh.push_back(&createHolder(0, 7.1, 10, 2, 2, 2.8, 0.4, N));
    mesh[8]->rotate(0, 1, 0, M_PI / 2);
    mesh[8]->setAmbientColor(AMB_OBSIDIAN);
    mesh[8]->setDiffuseColor(DIF_OBSIDIAN);
    mesh[8]->setSpecularColor(SPE_OBSIDIAN);
    mesh[8]->setShininess(100 * SHI_OBSIDIAN);
    mesh[6]->addDependency(*mesh[8]);

    // Pivots
    mesh.push_back(&createDuoCylinder(0, 0, 0, 1, 0.5, 12, 0.5, 0.3, N));
    mesh[9]->rotate(0, 0, 1, M_PI / 2);
    mesh[9]->translate(0, 8.9, 0);
    mesh[9]->translate(-1, 0, 0);
    mesh[9]->setAmbientColor(AMB_COPPER);
    mesh[9]->setDiffuseColor(DIF_COPPER);
    mesh[9]->setSpecularColor(SPE_COPPER);
    mesh[9]->setShininess(100 * SHI_COPPER);

    mesh.push_back(&createDuoCylinder(0, 0, 0, 1, 0.5, 12, 0.5, 0.3, N));
    mesh[10]->rotate(0, 0, 1, -M_PI / 2);
    mesh[10]->translate(0, 8.9, 0);
    mesh[10]->translate(1, 0, 0);
    mesh[10]->setAmbientColor(AMB_COPPER);
    mesh[10]->setDiffuseColor(DIF_COPPER);
    mesh[10]->setSpecularColor(SPE_COPPER);
    mesh[10]->setShininess(100 * SHI_COPPER);

    // Middle oval
    mesh.push_back(&createSpeacialOval(0, 8.1, 0, 0.5, 1, 10, 1.5, N));

    mesh[11]->setAmbientColor(AMB_COPPER);
    mesh[11]->setDiffuseColor(DIF_COPPER);
    mesh[11]->setSpecularColor(SPE_COPPER);
    mesh[11]->setShininess(100 * SHI_COPPER);
    mesh[11]->addDependency(*mesh[9]);
    mesh[11]->addDependency(*mesh[10]);
    mesh[3]->addDependency(*mesh[11]);

    // Latch
    mesh.push_back(&createDuoCylinder(0, 8.1, 4, 1.5, 0.2, 0, 0.5, 0.25, N));
    mesh[12]->setAmbientColor(AMB_RUBY);
    mesh[12]->setDiffuseColor(DIF_RUBY);
    mesh[12]->setSpecularColor(SPE_RUBY);
    mesh[12]->setShininess(100);
    mesh[4]->addDependency(*mesh[12]);

    // Initial position
    for (int i = 0; i < 100; i++)
        myKeyboard('5', 0, 0);
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
    glutCreateWindow("Vu Hoang Van - 1614063");
    init();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);

    glutKeyboardFunc(myKeyboard);
    glutSpecialFunc(mySpecialKeyboard);
    glutMainLoop();
    return 0;
}

#include "mesh.h"
#include <iomanip>
#define ZERO 0.000000001

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

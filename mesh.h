// #ifndef MESH_H
// #define MESH_H

// #include <vector>
// #include <cmath>
// #include <GL/glut.h>
// #include <string>
// #include <fstream>
// #include <iostream>

// #define ZERO 0.000000001

// using namespace std;

// class Vector;
// class Point;
// class Face;
// class Mesh;

// class Vector
// {

// public:
//     double x, y, z;
//     Vector()
//     {
//         x = y = z = 0;
//     }
//     Vector(double _x, double _y, double _z)
//     {
//         x = _x;
//         y = _y;
//         z = _z;
//     }

//     Vector(const Vector &vect)
//     {
//         x = vect.x;
//         y = vect.y;
//         z = vect.z;
//     }

//     Vector operator+(Vector &vect);
//     Vector operator-(Vector &vect);
//     Vector operator-() const;
//     Vector operator*(double factor);
//     Vector operator*(double factor) const;
//     Vector operator/(double divisor);
//     double operator*(Vector &vect);
//     Vector operator|(Vector &vect);
//     double len();
//     Vector normalize();
//     void print();
//     bool isZero();
// };

// class Point : public Vector
// {
// public:
//     Point() : Vector() {}

//     Point(double _x, double _y, double _z) : Vector(_x, _y, _z){};

//     Point(const Vector &&v)
//     {
//         x = v.x;
//         y = v.y;
//         z = v.z;
//     }
// };

// class Face
// {
// public:
//     vector<Point *> *mainVertice;
//     Vector norm;
//     vector<int> vertice;
//     Mesh *mesh;

//     float cs[4];
//     float cd[4];
//     float ca[4];
//     float shininess = 100.0;

//     Face()
//     {
//         mainVertice = NULL;
//         mesh = NULL;
//     }

//     Face(vector<Point *> *mainVertices) : Face()
//     {
//         mainVertice = mainVertices;
//     }

//     Face(const Face &face)
//     {
//         mainVertice = face.mainVertice;
//         norm = face.norm;
//         vertice = face.vertice;
//         mesh = face.mesh;
//         for (int i = 0; i < 4; i++)
//         {
//             cs[i] = face.cs[i];
//             cd[i] = face.cd[i];
//             ca[i] = face.ca[i];
//         }
//         shininess = face.shininess;
//     }
//     void setMainVertices(vector<Point *> &mainVertices);
//     void setAmbientColor(double red, double green, double blue, double alpha);
//     void setSpecularColor(double red, double green, double blue, double alpha);
//     void setDiffuseColor(double red, double green, double blue, double alpha);
//     void setShininess(double value);
//     void addPoint(int index);
//     void recalNorm();
//     void drawBorder();
//     void drawColor();
//     void setMesh(Mesh *m);
// };

// class Mesh
// {
// public:
//     bool reverse = false;
//     vector<Point *> mainVertice;
//     vector<Face *> face;
//     Mesh *dpO = NULL, *dpOX = NULL, *dpOY = NULL, *dpOZ = NULL;
//     vector<Mesh *> dp, dpX, dpY, dpZ;
//     Vector lightDirection;
//     double csr, csg, csb;
//     Mesh()
//     {
//         lightDirection = Vector(0, 0, 0);
//         csr = 1, csg = 1, csb = 1;
//     }
//     Mesh(string filename) : Mesh()
//     {
//         // Read data from a file

//         fstream f;
//         f.open(filename);
//         int nVertex = 0;
//         int nFace = 0;
//         f >> nVertex >> nFace;
//         mainVertice.clear();
//         for (int i = 0; i < nVertex; i++)
//         {
//             double _x, _y, _z;
//             f >> _x >> _y >> _z;
//             mainVertice.push_back(new Point(_x, _y, _z));
//         }
//         for (int i = 0; i < nFace; i++)
//         {
//             newFace();
//             Face &lface = getLastFace();
//             int vCount;
//             f >> vCount;
//             while (vCount--)
//             {
//                 int x;
//                 f >> x;
//                 lface.addPoint(x);
//             }
//             double r, g, b, a;
//             f >> r >> g >> b >> a;
//             lface.setAmbientColor(r, g, b, a);
//         }
//     }
//     void addPoint(Point &point);
//     void addFace(Face &_face);
//     void newFace();
//     Face &getLastFace();
//     void drawBorder();
//     void drawColor();
//     void translate(double dx, double dy, double dz);
//     void rotate(double ux, double uy, double uz, double angle);
//     void scale(double fx, double fy, double fz);
//     void addDependency(Mesh &m);
//     void writeToDisk(string filename);

//     void setAmbientColor(double red, double green, double blue, double alpha);
//     void setSpecularColor(double red, double green, double blue, double alpha);
//     void setDiffuseColor(double red, double green, double blue, double alpha);
//     void setShininess(double value);
// };

// Mesh &createHolder(double bottomX, double bottomY, double bottomZ,
//                    double lx, double lz, double height, double inR,
//                    int acc)
// {
//     Point *fdr = new Point(lx / 2, 0, -lz / 2);
//     Point *fdl = new Point(-lx / 2, 0, -lz / 2);
//     Point *bdl = new Point(-lx / 2, 0, lz / 2);
//     Point *bdr = new Point(lx / 2, 0, lz / 2);

//     Point *fur = new Point(lx / 2, height - lx / 2, -lz / 2);
//     Point *ful = new Point(-lx / 2, height - lx / 2, -lz / 2);
//     Point *bul = new Point(-lx / 2, height - lx / 2, lz / 2);
//     Point *bur = new Point(lx / 2, height - lx / 2, lz / 2);

//     Point *fir = new Point(inR, height - lx / 2, -lz / 2);
//     Point *fil = new Point(-inR, height - lx / 2, -lz / 2);
//     Point *bil = new Point(-inR, height - lx / 2, lz / 2);
//     Point *bir = new Point(inR, height - lx / 2, lz / 2);

//     Point *fjr = new Point(inR, 0, -lz / 2);
//     Point *fjl = new Point(-inR, 0, -lz / 2);
//     Point *bjl = new Point(-inR, 0, lz / 2);
//     Point *bjr = new Point(inR, 0, lz / 2);

//     int &N = acc;
//     map<void *, int> pMap;
//     Mesh &mesh = *(new Mesh());
//     int index = 0;

//     {
//         mesh.addPoint(*fdr);
//         pMap[fdr] = index++;

//         mesh.addPoint(*fdl);
//         pMap[fdl] = index++;

//         mesh.addPoint(*bdl);
//         pMap[bdl] = index++;

//         mesh.addPoint(*bdr);
//         pMap[bdr] = index++;

//         mesh.addPoint(*ful);
//         pMap[ful] = index++;

//         mesh.addPoint(*fur);
//         pMap[fur] = index++;

//         mesh.addPoint(*bul);
//         pMap[bul] = index++;

//         mesh.addPoint(*bur);
//         pMap[bur] = index++;

//         mesh.addPoint(*fil);
//         pMap[fil] = index++;

//         mesh.addPoint(*fir);
//         pMap[fir] = index++;

//         mesh.addPoint(*bil);
//         pMap[bil] = index++;

//         mesh.addPoint(*bir);
//         pMap[bir] = index++;

//         mesh.addPoint(*fjl);
//         pMap[fjl] = index++;

//         mesh.addPoint(*fjr);
//         pMap[fjr] = index++;

//         mesh.addPoint(*bjl);
//         pMap[bjl] = index++;

//         mesh.addPoint(*bjr);
//         pMap[bjr] = index++;
//     }

//     vector<Point *> fi, fo, bi, bo, fd, bd;
//     for (int i = 0; i < (int)N; i++)
//     {
//         double angle = 2.0 * M_PI / N * i;
//         Point *pf = new Point(inR * cos(angle), height - lx / 2 + inR * sin(angle), -lz / 2);
//         Point *pb = new Point(inR * cos(angle), height - lx / 2 + inR * sin(angle), lz / 2);
//         fi.push_back(pf);
//         bi.push_back(pb);
//         mesh.addPoint(*pf);
//         mesh.addPoint(*pb);
//         pMap[pf] = index++;
//         pMap[pb] = index++;
//     }
//     fi.push_back(fi[0]);
//     bi.push_back(bi[0]);

//     for (int i = N / 2; i < (int)N + 1; i++)
//     {
//         double angle = 2.0 * M_PI / N * i;
//         Point *pf = new Point(inR * cos(angle), 0, -lz / 2);
//         Point *pb = new Point(inR * cos(angle), 0, lz / 2);
//         fd.push_back(pf);
//         bd.push_back(pb);
//         mesh.addPoint(*pf);
//         mesh.addPoint(*pb);
//         pMap[pf] = index++;
//         pMap[pb] = index++;
//     }

//     double outR = lx / 2;
//     for (int i = 0; i < (int)N / 2 + 1; i++)
//     {
//         double angle = 2.0 * M_PI / N * i;
//         Point *pf = new Point(outR * cos(angle), height - lx / 2 + outR * sin(angle), -lz / 2);
//         Point *pb = new Point(outR * cos(angle), height - lx / 2 + outR * sin(angle), lz / 2);
//         fo.push_back(pf);
//         bo.push_back(pb);
//         mesh.addPoint(*pf);
//         mesh.addPoint(*pb);
//         pMap[pf] = index++;
//         pMap[pb] = index++;
//     }
//     fo.push_back(fo[0]);
//     bo.push_back(bo[0]);

//     // Bottom face
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[fdr]);
//         face.addPoint(pMap[bdr]);
//         face.addPoint(pMap[bdl]);
//         face.addPoint(pMap[fdl]);
//     }

//     // Left face
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[fdl]);
//         face.addPoint(pMap[bdl]);
//         face.addPoint(pMap[bul]);
//         face.addPoint(pMap[ful]);
//     }

//     // Right face
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[fdr]);
//         face.addPoint(pMap[fur]);
//         face.addPoint(pMap[bur]);
//         face.addPoint(pMap[bdr]);
//     }

//     // Front face 1
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[fdl]);
//         face.addPoint(pMap[ful]);
//         face.addPoint(pMap[fil]);
//         face.addPoint(pMap[fjl]);
//     }

//     // Front face 2
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[fdr]);
//         face.addPoint(pMap[fjr]);
//         face.addPoint(pMap[fir]);
//         face.addPoint(pMap[fur]);
//     }

//     // Back face 1
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[bdl]);
//         face.addPoint(pMap[bjl]);
//         face.addPoint(pMap[bil]);
//         face.addPoint(pMap[bul]);
//     }

//     // Back face 2
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[bdr]);
//         face.addPoint(pMap[bur]);
//         face.addPoint(pMap[bir]);
//         face.addPoint(pMap[bjr]);
//     }
//     // Outside face
//     for (int i = 0; i < (int)N / 2; i++)
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[fo[i]]);
//         face.addPoint(pMap[fo[i + 1]]);
//         face.addPoint(pMap[bo[i + 1]]);
//         face.addPoint(pMap[bo[i]]);
//     }

//     // Inside face
//     for (int i = 0; i < (int)N; i++)
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[fi[i]]);
//         face.addPoint(pMap[bi[i]]);
//         face.addPoint(pMap[bi[(i + 1) % N]]);
//         face.addPoint(pMap[fi[(i + 1) % N]]);
//     }

//     // Front face
//     for (int i = 0; i < (int)N / 2; i++)
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[fi[i]]);
//         face.addPoint(pMap[fi[i + 1]]);
//         face.addPoint(pMap[fo[i + 1]]);
//         face.addPoint(pMap[fo[i]]);
//     }

//     for (int i = N / 2; i < (int)N; i++)
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[fi[i]]);
//         face.addPoint(pMap[fi[i + 1]]);
//         face.addPoint(pMap[fd[i - N / 2 + 1]]);
//         face.addPoint(pMap[fd[i - N / 2]]);
//     }

//     // Back face
//     for (int i = 0; i < (int)N / 2; i++)
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[bi[i]]);
//         face.addPoint(pMap[bo[i]]);
//         face.addPoint(pMap[bo[i + 1]]);
//         face.addPoint(pMap[bi[i + 1]]);
//     }

//     for (int i = N / 2; i < (int)N - 1; i++)
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[bi[i]]);
//         face.addPoint(pMap[bd[i - N / 2]]);
//         face.addPoint(pMap[bd[i - N / 2 + 1]]);
//         face.addPoint(pMap[bi[i + 1]]);
//     }
//     mesh.translate(bottomX, bottomY, bottomZ);

//     // Clean up
//     delete fdr;
//     delete fdl;
//     delete bdl;
//     delete bdr;

//     delete fur;
//     delete ful;
//     delete bul;
//     delete bur;

//     delete fir;
//     delete fil;
//     delete bil;
//     delete bir;

//     delete fjr;
//     delete fjl;
//     delete bjl;
//     delete bjr;

//     pMap.clear();
//     for (int i = 0; i < (int)fi.size() - 1; i++)
//         delete fi[i];
//     for (int i = 0; i < (int)fo.size() - 1; i++)
//         delete fo[i];
//     for (int i = 0; i < (int)fd.size() - 1; i++)
//         delete fd[i];
//     for (int i = 0; i < (int)bi.size() - 1; i++)
//         delete bi[i];
//     for (int i = 0; i < (int)bo.size() - 1; i++)
//         delete bo[i];
//     for (int i = 0; i < (int)bd.size() - 1; i++)
//         delete bd[i];
//     fi.clear();
//     fo.clear();
//     fd.clear();
//     bi.clear();
//     bo.clear();
//     bd.clear();

//     return mesh;
// }

// Mesh &createCylinder(double bottomX, double bottomY, double bottomZ,
//                      double radius, double height,
//                      int acc)
// {

//     int &N = acc;
//     Mesh &mesh = *(new Mesh("empty.dat"));
//     Point *topCenter = new Point(0, height, 0);
//     Point *bottomCenter = new Point(0, 0, 0);
//     mesh.addPoint(*topCenter);
//     mesh.addPoint(*bottomCenter);

//     map<void *, int> pMap;
//     pMap[topCenter] = 0;
//     pMap[bottomCenter] = 1;
//     int index = 2;
//     Point *topPoint = new Point[N]();
//     Point *bottomPoint = new Point[N]();

//     for (int i = 0; i < N; i++)
//     {
//         double angle = 2.0 * M_PI / N * i;
//         topPoint[i] = Point(radius * cos(angle), height, radius * sin(angle));
//         bottomPoint[i] = Point(radius * cos(angle), 0, radius * sin(angle));
//         mesh.addPoint(topPoint[i]);
//         mesh.addPoint(bottomPoint[i]);
//         pMap[topPoint + i] = index++;
//         pMap[bottomPoint + i] = index++;
//     }

//     for (int i = N - 1; i >= 0; i--)
//     {
//         // Top face
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[topPoint + (i + 1) % N]);
//         face.addPoint(pMap[topPoint + i]);

//         face.addPoint(0);
//     }

//     for (int i = 0; i < N; i++)
//     {

//         {
//             // Bottom face
//             mesh.newFace();
//             Face &face = mesh.getLastFace();
//             face.addPoint(pMap[bottomPoint + i]);
//             face.addPoint(pMap[bottomPoint + (i + 1) % N]);
//             face.addPoint(1);
//         }
//         {
//             // Side face
//             mesh.newFace();
//             Face &face = mesh.getLastFace();
//             face.addPoint(pMap[topPoint + i]);
//             face.addPoint(pMap[topPoint + (i + 1) % N]);
//             face.addPoint(pMap[bottomPoint + (i + 1) % N]);
//             face.addPoint(pMap[bottomPoint + i]);
//         }
//     }

//     mesh.translate(bottomX, bottomY, bottomZ);

//     // Clean up
//     delete[] topPoint;
//     delete[] bottomPoint;
//     delete topCenter;
//     delete bottomCenter;

//     return mesh;
// }

// Mesh &createRectangular(double bottomX, double bottomY, double bottomZ,
//                         double rx, double ry, double rz)
// {
//     Point *bsw = new Point(-rx / 2, 0, -rz / 2);
//     Point *bse = new Point(rx / 2, 0, -rz / 2);
//     Point *bnw = new Point(-rx / 2, 0, rz / 2);
//     Point *bne = new Point(rx / 2, 0, rz / 2);

//     Point *usw = new Point(-rx / 2, ry, -rz / 2);
//     Point *use = new Point(rx / 2, ry, -rz / 2);
//     Point *unw = new Point(-rx / 2, ry, rz / 2);
//     Point *une = new Point(rx / 2, ry, rz / 2);

//     map<void *, int> pMap;
//     Mesh &mesh = *(new Mesh());
//     int index = 0;

//     mesh.addPoint(*bsw);
//     pMap[bsw] = index++;

//     mesh.addPoint(*bse);
//     pMap[bse] = index++;

//     mesh.addPoint(*bnw);
//     pMap[bnw] = index++;

//     mesh.addPoint(*bne);
//     pMap[bne] = index++;

//     mesh.addPoint(*usw);
//     pMap[usw] = index++;

//     mesh.addPoint(*use);
//     pMap[use] = index++;

//     mesh.addPoint(*unw);
//     pMap[unw] = index++;

//     mesh.addPoint(*une);
//     pMap[une] = index++;
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[usw]);
//         face.addPoint(pMap[use]);
//         face.addPoint(pMap[bse]);
//         face.addPoint(pMap[bsw]);
//     }
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[use]);
//         face.addPoint(pMap[une]);
//         face.addPoint(pMap[bne]);
//         face.addPoint(pMap[bse]);
//     }
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[une]);
//         face.addPoint(pMap[unw]);
//         face.addPoint(pMap[bnw]);
//         face.addPoint(pMap[bne]);
//     }
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[unw]);
//         face.addPoint(pMap[usw]);
//         face.addPoint(pMap[bsw]);
//         face.addPoint(pMap[bnw]);
//     }
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[unw]);
//         face.addPoint(pMap[une]);
//         face.addPoint(pMap[use]);
//         face.addPoint(pMap[usw]);
//     }
//     {
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[bnw]);
//         face.addPoint(pMap[bsw]);
//         face.addPoint(pMap[bse]);
//         face.addPoint(pMap[bne]);
//     }
//     mesh.translate(bottomX, bottomY, bottomZ);
//     return mesh;
// }

// Mesh &createDuoCylinder(double bottomX, double bottomY, double bottomZ,
//                         double h1, double hs, double h2, double rd, double ru,
//                         int acc)
// {
//     int &N = acc;
//     Mesh &mesh = *(new Mesh("empty.dat"));
//     Point *bottomCenter = new Point(0, 0, 0);
//     Point *topCenter = new Point(0, h1 + hs + h2, 0);
//     mesh.addPoint(*bottomCenter);
//     mesh.addPoint(*topCenter);
//     map<void *, int> pMap;
//     pMap[bottomCenter] = 0;
//     pMap[topCenter] = 1;
//     int index = 2;
//     Point *l0 = new Point[N]();
//     Point *l1 = new Point[N]();
//     Point *l2 = new Point[N]();
//     Point *l3 = new Point[N]();
//     for (int i = 0; i < N; i++)
//     {
//         double angle = 2.0 * M_PI / N * i;
//         l0[i] = Point(rd * cos(angle), 0, rd * sin(angle));
//         l1[i] = Point(rd * cos(angle), h1, rd * sin(angle));
//         l2[i] = Point(ru * cos(angle), h1 + hs, ru * sin(angle));
//         l3[i] = Point(ru * cos(angle), h1 + hs + h2, ru * sin(angle));
//         mesh.addPoint(l0[i]);
//         mesh.addPoint(l1[i]);
//         mesh.addPoint(l2[i]);
//         mesh.addPoint(l3[i]);
//         pMap[l0 + i] = index++;
//         pMap[l1 + i] = index++;
//         pMap[l2 + i] = index++;
//         pMap[l3 + i] = index++;
//     }

//     for (int i = 0; i < N; i++)
//     {
//         {
//             // Top face
//             mesh.newFace();
//             Face &face = mesh.getLastFace();
//             face.addPoint(pMap[l3 + i]);
//             face.addPoint(pMap[topCenter]);
//             face.addPoint(pMap[l3 + (i + 1) % N]);
//         }
//         {
//             // Bottom face
//             mesh.newFace();
//             Face &face = mesh.getLastFace();
//             face.addPoint(pMap[l0 + (i + 1) % N]);
//             face.addPoint(pMap[bottomCenter]);
//             face.addPoint(pMap[l0 + i]);
//         }
//         {
//             // Bottom side face
//             mesh.newFace();
//             Face &face = mesh.getLastFace();
//             face.addPoint(pMap[l1 + i]);
//             face.addPoint(pMap[l1 + (i + 1) % N]);
//             face.addPoint(pMap[l0 + (i + 1) % N]);
//             face.addPoint(pMap[l0 + i]);
//         }
//         {
//             // Middle side face
//             mesh.newFace();
//             Face &face = mesh.getLastFace();
//             face.addPoint(pMap[l2 + i]);
//             face.addPoint(pMap[l2 + (i + 1) % N]);
//             face.addPoint(pMap[l1 + (i + 1) % N]);
//             face.addPoint(pMap[l1 + i]);
//         }
//         {
//             // Top side face
//             mesh.newFace();
//             Face &face = mesh.getLastFace();
//             face.addPoint(pMap[l3 + i]);
//             face.addPoint(pMap[l3 + (i + 1) % N]);
//             face.addPoint(pMap[l2 + (i + 1) % N]);
//             face.addPoint(pMap[l2 + i]);
//         }
//     }

//     mesh.translate(bottomX, bottomY, bottomZ);
//     // Clean up
//     delete[] l0;
//     delete[] l1;
//     delete[] l2;
//     delete[] l3;
//     // delete topCenter;
//     // delete bottomCenter;
//     return mesh;
// }

// Mesh &createSpeacialOval(double bottomX, double bottomY, double bottomZ,
//                          double ri, double ro, double len, double height,
//                          int acc)
// {
//     int &N = acc;
//     Mesh &mesh = *(new Mesh("empty.dat"));
//     map<void *, int> pMap;
//     Point *_ui = new Point[N + 2]();
//     Point *_uo = new Point[N + 2]();
//     Point *_di = new Point[N + 2]();
//     Point *_do = new Point[N + 2]();
//     int index = 0;
//     for (int i = 0; i <= N / 2; i++)
//     {
//         double angle = 2.0 * M_PI / N * i;
//         _ui[i] = Point(ri * cos(angle), height, ri * sin(angle) + len / 2);
//         mesh.addPoint(_ui[i]);
//         pMap[_ui + i] = index++;
//         _uo[i] = Point(ro * cos(angle), height, ro * sin(angle) + len / 2);
//         mesh.addPoint(_uo[i]);
//         pMap[_uo + i] = index++;
//         _di[i] = Point(ri * cos(angle), 0, ri * sin(angle) + len / 2);
//         mesh.addPoint(_di[i]);
//         pMap[_di + i] = index++;
//         _do[i] = Point(ro * cos(angle), 0, ro * sin(angle) + len / 2);
//         mesh.addPoint(_do[i]);
//         pMap[_do + i] = index++;
//     }
//     for (int i = N / 2 + 1; i <= N + 1; i++)
//     {
//         double angle = 2.0 * M_PI / N * (i - 1);
//         _ui[i] = Point(ri * cos(angle), height, ri * sin(angle) - len / 2);
//         mesh.addPoint(_ui[i]);
//         pMap[_ui + i] = index++;
//         _uo[i] = Point(ro * cos(angle), height, ro * sin(angle) - len / 2);
//         mesh.addPoint(_uo[i]);
//         pMap[_uo + i] = index++;
//         _di[i] = Point(ri * cos(angle), 0, ri * sin(angle) - len / 2);
//         mesh.addPoint(_di[i]);
//         pMap[_di + i] = index++;
//         _do[i] = Point(ro * cos(angle), 0, ro * sin(angle) - len / 2);
//         mesh.addPoint(_do[i]);
//         pMap[_do + i] = index++;
//     }

//     for (int i = 0; i < N + 2; i++)
//     {
//         // Top face
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[_ui + i]);
//         face.addPoint(pMap[_ui + (i + 1) % (N + 2)]);
//         face.addPoint(pMap[_uo + (i + 1) % (N + 2)]);
//         face.addPoint(pMap[_uo + i]);
//     }
//     for (int i = 0; i < N + 2; i++)
//     {
//         // Bottom face
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[_di + i]);
//         face.addPoint(pMap[_do + i]);
//         face.addPoint(pMap[_do + (i + 1) % (N + 2)]);
//         face.addPoint(pMap[_di + (i + 1) % (N + 2)]);
//     }
//     for (int i = 0; i < N + 2; i++)
//     {
//         // Inside face
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[_ui + i]);
//         face.addPoint(pMap[_di + i]);
//         face.addPoint(pMap[_di + (i + 1) % (N + 2)]);
//         face.addPoint(pMap[_ui + (i + 1) % (N + 2)]);
//     }
//     for (int i = 0; i < N + 2; i++)
//     {
//         // Outside face
//         mesh.newFace();
//         Face &face = mesh.getLastFace();
//         face.addPoint(pMap[_uo + i]);
//         face.addPoint(pMap[_uo + (i + 1) % (N + 2)]);
//         face.addPoint(pMap[_do + (i + 1) % (N + 2)]);
//         face.addPoint(pMap[_do + i]);
//     }

//     mesh.translate(bottomX, bottomY, bottomZ);
//     delete[] _uo;
//     delete[] _ui;
//     delete[] _do;
//     delete[] _di;
//     return mesh;
// }

// #endif
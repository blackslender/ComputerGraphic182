// #include "mesh.h"
// #include <map>

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
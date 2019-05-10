// #include "mesh.h"
// #include <map>
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
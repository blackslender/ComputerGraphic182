// #include "mesh.h"
// #include <map>
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
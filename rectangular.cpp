// #include "mesh.h"
// #include <map>

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
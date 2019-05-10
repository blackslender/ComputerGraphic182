// #include <map>
// #include "mesh.h"
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
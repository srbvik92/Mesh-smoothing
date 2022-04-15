#include "io.h"
#include <iostream>
#include <fstream>
#include <cfloat>
#include <string>
#include "LinearSolver/cghs.h"
#include "LinearSolver/cblas.h"

#ifdef __APPLE__
#define MAX std::numeric_limits<double>::max()
#else
//#include <values.h>
#define MAX DBL_MAX
#endif

struct Operator {
    int n;
    double lambda, dt;
    std::vector<Tri> tris;
};

std::vector<Eigen::Vector3d> laplacian(const std::vector<Eigen::Vector3d>& pts, const std::vector<Tri>& tris);

int existvert(Eigen::Vector3d pt, std::vector<Eigen::Vector3d> points);
bool existtri(Tri t, std::vector<Tri> tris, std::vector<Eigen::Vector3d> points);

//this function based on: https://www.csee.umbc.edu/~adamb/435/implicit.cpp
void mult(const Operator& op, double* v, double* w) {
    std::vector<Eigen::Vector3d> pts;
    pts.resize(op.n);
    for (unsigned int i = 0; i < pts.size(); i++) {
        //printf("flag1 i=%d\n", i);
        pts[i] = Eigen::Vector3d(v[3 * i + 0], v[3 * i + 1], v[3 * i + 2]);
    }

    std::vector<Eigen::Vector3d> l = laplacian(pts, op.tris);
    for (unsigned int i = 0; i < pts.size(); i++) {
        l[i] *= op.lambda * op.dt;
        w[3 * i + 0] = v[3 * i + 0] - l[i][0];
        w[3 * i + 1] = v[3 * i + 1] - l[i][1];
        w[3 * i + 2] = v[3 * i + 2] - l[i][2];
    }

}

//laplacian
std::vector<Eigen::Vector3d> laplacian(const std::vector<Eigen::Vector3d>& points, const std::vector<Tri>& tris) {
    std::vector<Eigen::Vector3d> laplac;

    // initialize accumulators
    for (unsigned int i = 0; i < points.size(); i++) {
        laplac.push_back(Eigen::Vector3d(0.0, 0.0, 0.0));
    }

    // loop over triangles
    for (unsigned int i = 0; i < tris.size(); i++) {
        const Tri& t = tris[i];
        for (int c = 0; c < 3; c++) {
            laplac[t.indices[c]] += points[t.indices[(c + 1) % 3]] + points[t.indices[(c + 2) % 3]] - 2 * points[t.indices[c]];
        }
    }
    return laplac;
}


bool readObjFile(const char *fname, std::vector<Eigen::Vector3d> &pts, std::vector<Tri> &triangles)
{
  char c[500];
  Eigen::Vector3d lc(MAX,MAX,MAX), uc(-MAX,-MAX,-MAX);

  int numVertices=0, numFaces=0;
  bool normals = false, texture = false;
  int tint;
  char ch;
  int p, q, r;
  double x, y, z;
  std::vector<Eigen::Vector3d>::iterator v;
  std::vector<Tri>::iterator t;
  std::ifstream in1(fname, std::ios::in);
  if (!in1.is_open()) {
    return false;
  }
  in1.flags(in1.flags() & ~std::ios::skipws);

  while (in1>>ch) {
    if (ch == 'v') {
      in1>>ch;
      if (ch == ' ') numVertices++;
      else if (ch == 'n') normals = true;
      else if (ch == 't') texture = true;
      else std::cerr<<"error \'"<<ch<<"\'"<<std::endl;
    } else if (ch == '#') {
      while (in1 >> ch && ch != '\n') ; // Read to the end of the line.
    } else if (ch == 'f') numFaces++;
  }
  in1.close();

  pts.resize(numVertices);
  triangles.resize(numFaces);
  v = pts.begin();
  t = triangles.begin();

  std::ifstream in(fname, std::ios::in);
  if (!in.is_open()) {
    return false;
  }

  while (in>>ch) {
    if (ch == '#') {
      in.getline(c,500);
      continue;
    }
    if (ch == 'g') {
      in.getline(c,500);
      continue;
    }
    if (ch == 's') {
      in.getline(c,500);
      continue;
    }
    if (ch == 'm') {
      in.getline(c,500);
      continue;
    }
    if (ch == 'u') {
      in.getline(c,500);
      continue;
    }
    if (ch == 'v') {
      ch = in.peek();
      if (ch != 't' && ch != 'n') {
	in>>x>>y>>z;
	(*v)<<x,y,z;
	if ((*v)[0] < lc[0]) lc[0] = (*v)[0];
	if ((*v)[1] < lc[1]) lc[1] = (*v)[1];
	if ((*v)[2] < lc[2]) lc[2] = (*v)[2];
	if ((*v)[0] > uc[0]) uc[0] = (*v)[0];
	if ((*v)[1] > uc[1]) uc[1] = (*v)[1];
	if ((*v)[2] > uc[2]) uc[2] = (*v)[2];
	v++;
      } else {
	in.getline(c, 500);
      }
      continue;
    }
    if (ch == 'f') {
      if (normals && texture) {
	in>>p>>ch>>tint>>ch>>tint>>q>>ch>>tint>>ch>>tint>>r>>ch>>tint>>ch>>tint;
      } else if (normals) {
	in>>p>>ch>>ch>>tint>>q>>ch>>ch>>tint>>r>>ch>>ch>>tint;
      } else if (texture) {
	in>>p>>ch>>tint>>q>>ch>>tint>>r>>ch>>tint;
      } else {
	in>>p>>q>>r;
      }
      (*t)[0] = p-1;
      (*t)[1] = q-1;
      (*t)[2] = r-1;
      t++;
      continue;
    }
  }
  in.close();
  return true;
}


void writeObjFile(const char *fname, const std::vector<Eigen::Vector3d> &meshPts, const std::vector<Tri> &triangles) {
  std::ofstream out;
  std::vector<Eigen::Vector3d>::const_iterator p;
  std::vector<Tri>::const_iterator t;

  out.open(fname);

  for (p=meshPts.begin(); p!=meshPts.end(); p++) 
    out<<"v "<<(*p)[0]<<" "<<(*p)[1]<<" "<<(*p)[2]<<std::endl;
  
  for (t=triangles.begin(); t!=triangles.end(); t++) 
    out<<"f "<<(*t)[0]+1<<" "<<(*t)[1]+1<<" "<<(*t)[2]+1<<std::endl;

  out.close();
}


int main(int argc, char* argv[]) {

  //basic variables
  char* inputfile, * outputfile;
  double lambda, dt;
  int iterations, subdpass;
  double eps = 0;

  // read command line args
  if (argc < 2) {
      //usage();
      std::cout << "too few arguments..... exiting";
      exit(1);
  }

  bool cotangent = std::string((argv[1])) == "-c";        //for cotangent weights
  bool gvc = std::string((argv[1])) == "-v";              //for global volume preservation
  bool biharmonic = std::string((argv[1])) == "-b";       //for biharmonic operator
  bool loopsubd = std::string((argv[1])) == "-s";         //for loop subdivision
  bool implicit = std::string((argv[1])) == "-i";

  //std::cout << cotangent << gvc << biharmonic << loopsubd << implicit;

  //if not implicit
  /*if (!implicit && argc != 6) {

    //if not loopsubd
    if(!loopsubd && argc != 6){
      std::cout << "Exactly 5 args expected!" << std::endl;
      exit(1);
    }
    
  } */

  /*else if(!loopsubd && argc != 6) {
    std::cout << "Exactly 5 args expected!" << std::endl;
    exit(1);
  } */

  //for normal smoothing i.e. default case'
  /*if(!cotangent && !implicit && !gvc && !biharmonic && !loopsubd && argc != 6){
    std::cout << "Exactly 5 args expected!" << std::endl;
    exit(1);
  } */

  /*if (!loopsubd && argc != 6) {
    std::cout << "Exactly 5 args expected!" << std::endl;
    exit(1);
  } */

  //read eps for implicit and volume preservation
  if (implicit || gvc) {
      eps = std::stod(argv[2]);
      //just advance so that we can parse the arguments the same way
      argv = &(argv[2]);
  }

  //read subdivision pass for loop subdivision
  if(loopsubd){
      subdpass = std::stoi(argv[2]);
      argv = &(argv[2]);
  }
  

  inputfile = argv[1];
  outputfile = argv[2];

  sscanf(argv[3], "%lf", &lambda);
  sscanf(argv[4], "%lf", &dt);
  sscanf(argv[5], "%d", &iterations);

  //two variables for points of meshe, one will keep original values and other will be modified
  std::vector<Eigen::Vector3d> points, newpoints;

  //for storing triangles of mesh
  std::vector<Tri> tris;

  //Laplacian amd m
  std::vector<float> m;
  std::vector<Eigen::Vector3d> L;

  readObjFile(inputfile, points, tris);   //reading the obj file;

  newpoints = points;

  //normal laplacian smoothing
  if(!cotangent && !implicit && !gvc && !biharmonic && !loopsubd) {

    std::cout << "Smoothing using laplacian\n";
    for (int i = 0; i < points.size(); i++) {
        m.push_back(0.0);
        Eigen::Vector3d a(0.0, 0.0, 0.0);
        L.push_back(a);
    }

    //iterating for number of triangles
    for (int itr = 0; itr < iterations; itr++) {
        
        //iterating for each triangle and calculating laplacian for each vertices of triangle
        for (int i = 0; i < tris.size(); i++) {
            Tri& tri = tris[i];

            //for 0th vertex of triangle
            L[tri[0]] += newpoints[tri[1]] - newpoints[tri[0]]
                + newpoints[tri[2]] - newpoints[tri[0]];

            m[tri[0]] += 2;

            //for 1st vertex of triangle
            L[tri[1]] += newpoints[tri[0]] - newpoints[tri[1]]
                + newpoints[tri[2]] - newpoints[tri[1]];

            m[tri[1]] += 2;

            //for 2nd vertex of triangle
            L[tri[2]] += newpoints[tri[1]] - newpoints[tri[2]]
                + newpoints[tri[0]] - newpoints[tri[2]];

            m[tri[2]] += 2;

        }

        //updating values of points based on laplacian calculation
        for (int i = 0; i < points.size(); i++) {
            newpoints[i] += lambda * dt * (L[i] / m[i]);
        }
   
    }
    writeObjFile(outputfile, newpoints, tris);                          //write the obj file for laplacian smoothing
    std::cout << "New obj file named- " << outputfile << " created";
  }


  //implicit integration based on professor's code
  if (implicit) {
    std::cout << "Implicit integration smoothing\n";

    newpoints = points;
    
    struct Operator op;
    op.n = newpoints.size(); op.lambda = lambda; op.dt = dt;
    op.tris = tris;

    //const int size = 3 * pts.size();

    double* b = new double[3 * newpoints.size()];
    double* x = new double[3 * newpoints.size()];
 
    for (unsigned int i = 0; i < newpoints.size(); i++) {
        for (int c = 0; c < 3; c++) {
            b[3 * i + c] = newpoints[i][c];
            x[3 * i + c] = newpoints[i][c];
        }
    }

    cghs<Operator>(3 * newpoints.size(), op, b, x, eps, true);

    for (unsigned int i = 0; i < newpoints.size(); i++) {
      newpoints[i] = Eigen::Vector3d(x[3 * i + 0], x[3 * i + 1], x[3 * i + 2]);
    }
    delete[] b;
    delete[] x;

    writeObjFile(outputfile, newpoints, op.tris);
    std::cout << "New obj file named- " << outputfile << " created";
  }

  //for volume computation
  if (gvc) {
    std::cout << "Global volume preservation smoothing\n";
    //calculate integration for volume first
    newpoints = points;
    
    struct Operator op;
    op.n = newpoints.size(); op.lambda = lambda; op.dt = dt;
    op.tris = tris;

    //const int size = 3 * pts.size();

    double* b = new double[3 * newpoints.size()];
    double* x = new double[3 * newpoints.size()];
 
    for (unsigned int i = 0; i < newpoints.size(); i++) {
        for (int c = 0; c < 3; c++) {
            b[3 * i + c] = newpoints[i][c];
            x[3 * i + c] = newpoints[i][c];
        }
    }

    cghs<Operator>(3 * newpoints.size(), op, b, x, eps, true);

    for (unsigned int i = 0; i < newpoints.size(); i++) {
        newpoints[i] = Eigen::Vector3d(x[3 * i + 0], x[3 * i + 1], x[3 * i + 2]);
    }
    delete[] b;
    delete[] x;

    //initial and volume after integration
    double V0, Vn  = 0;
    //newpoints = points;

    //calculating default volume
    for (int i = 0; i < tris.size(); i++) {
        Tri& tri = tris[i];
        
        Eigen::Vector3d gk = (points[tri[0]] + points[tri[1]] + points[tri[2]]);
        Eigen::Vector3d nk = ((points[tri[1]] - points[tri[0]]).cross(points[tri[2]] - points[tri[0]]));

        double temp = gk.dot(nk);
        V0 += temp;
    }
    V0 = V0/6.0;

    //calculating volume after integration
    for (int i = 0; i < tris.size(); i++) {
        Tri& tri = tris[i];

        Eigen::Vector3d gk = (newpoints[tri[0]] + newpoints[tri[1]] + newpoints[tri[2]]);
        Eigen::Vector3d nk = ((newpoints[tri[1]] - newpoints[tri[0]]).cross(newpoints[tri[2]] - newpoints[tri[0]]));

        double vn = gk.dot(nk);
        Vn = Vn + vn;
    }
    Vn = Vn / 6;

    double beta = pow((V0 / Vn), 1 / 3);            //calculating beta

    //modifying vertices after volume calculation by multiplying them with beta
    for (int i = 0; i < points.size(); i++) {
        newpoints[i] = beta * points[i];
    }
    
    writeObjFile(outputfile, newpoints, op.tris);                     //write obj file for volume preservation
    std::cout << "New obj file named- " << outputfile << " created\n";
  }

  //loop subdivision
  if(loopsubd) {
    
    std::cout << "Subdividing mesh \n";
    std::vector<Eigen::Vector3d> subdpoints;      //define new points and triangle variables
    std::vector<Tri> subdtris;

    const char* subd = "outsubd.obj";

    newpoints = points;

    for (int i = 0; i < tris.size(); i++) {
      Tri& tri = tris[i];
      int posx, posy, posz;

        //creating new vertices by taking average of points
        //for 0th and 1st vertex of triangle
        Eigen::Vector3d x = (newpoints[tri[0]] + newpoints[tri[1]]) / 2.0;

        //check whether vertex already exists in the list and if exists, return the position
        posx = existvert(x, subdpoints);            
        if (posx==-1) {
            subdpoints.push_back(x);
            posy = subdpoints.size();     //save the position of inserted vertex
        }
        
        //for 1st and 2ndvertex of triangle
        Eigen::Vector3d y = (newpoints[tri[1]] + newpoints[tri[2]]) / 2.0;

        //check whether vertex already exists in the list  and if exists, return the position
        posy = existvert(y, subdpoints);
        if (posy==-1) {
            subdpoints.push_back(y);
            posy = subdpoints.size();     //save the position of inserted vertex
        }

        //for 2nd and 0th vertex of triangle
        Eigen::Vector3d z = (newpoints[tri[0]] + newpoints[tri[2]]) / 2.0;

        //check whether vertex already exists in the list and if exists, return the position
        posz = existvert(z, subdpoints);
        if (posz==-1) {
            subdpoints.push_back(z);
            posz = subdpoints.size();     //save the position of inserted vertex
        }

        //create four triangles with possible combinations

        //0th, x and z vertices
        Tri tri0(tri[0], posx, posz);
        bool testtri = false;

        //checking whether triangle already exists in the list or not
        testtri = existtri(tri0, subdtris, subdpoints);

        //if doesn't exist then push it into list
        if (!testtri) {
          subdtris.push_back(tri0);
        }

        //1st, x and y vertices
        Tri tri1(tri[1], posx, posy);
        bool testtri1 = false;

        //checking whether triangle already exists in the list or not
        testtri1 = existtri(tri1, subdtris, subdpoints);

        //if doesn't exist then push it into list
        if (!testtri1) {
            subdtris.push_back(tri1);
        }

        //2nd, y and z vertices
        Tri tri2(tri[2], posx, posy);
        bool testtri2 = false;

        //checking whether triangle already exists in the list or not
        testtri2 = existtri(tri2, subdtris, subdpoints);

        //if doesn't exist then push it into list
        if (!testtri2) {
            subdtris.push_back(tri2);
        }

        //x, y and z vertices
        Tri tri3(posx, posy, posz);
        bool testtri3 = false;

        //checking whether triangle already exists in the list or not
        testtri3 = existtri(tri3, subdtris, subdpoints);

        //if doesn't exist then push it into list
        if (!testtri3) {
            subdtris.push_back(tri3);
        }

    }

    writeObjFile(outputfile, subdpoints, subdtris);                   //write the output subdivided file
    std::cout << "New obj file named- " << outputfile << " created";

  }
  return 1;
} 


//check if vertices exist in the list or not
int existvert(Eigen::Vector3d pt, std::vector<Eigen::Vector3d> points) {
  //if the vertex point is found in the list then return the position or return -1
  for (int i = 0; i < points.size(); i++) {
    if (pt == points[i]) {
      return i;
    }
  }
  return -1;
}

//check if triangle exist in the list or not
bool existtri(Tri t, std::vector<Tri> tris, std::vector<Eigen::Vector3d> points) {
  for (int i = 0; i < tris.size(); i++) {
    Tri& tri = tris[i];                                         

      //if all the vertices of triangle matches with any other triangle, then return true, else return false
      if (t[0] == tri[0] && t[1] == tri[1] && t[2] == tri[2]) {
        return true;
      }
    }
  return false;
}
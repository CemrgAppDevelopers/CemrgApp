/*=========================================================================

Program:   Medical Imaging & Interaction Toolkit
Language:  C++
Date:      $Date$
Version:   $Revision$

Copyright (c) German Cancer Research Center, Division of Medical and
Biological Informatics. All rights reserved.
See MITKCopyright.txt or http://www.mitk.org/copyright.html for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
/*=========================================================================
CEMRG CMD CARP UTILS
This app serves as a template for the command line apps to be implemented
in the framework.
=========================================================================*/

// Qmitk
#include <mitkIOUtil.h>
#include <mitkSurface.h>
#include <mitkImage.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkCommandLineParser.h>
#include <mitkImagePixelReadAccessor.h>

// VTK
#include <vtkImplicitBoolean.h>
#include <vtkImplicitVolume.h>
#include <vtkImplicitDataSet.h>
#include <vtkClipPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkPlane.h>
#include <vtkSphere.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkMath.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPolyDataNormals.h>
#include <vtkIdList.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkDataSetSurfaceFilter.h>

// ITK
#include <itkPoint.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkOrientImageFilter.h>

// Qt
#include <QtDebug>
#include <QDir>
#include <QString>
#include <QFileInfo>
#include <QProcess>
#include <QMessageBox>
#include <numeric>

#include <CemrgScar3D.h>
#include <CemrgCommandLine.h>

#include <algorithm>
#include <string>
// #include <math.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

typedef struct {
  double s;
  double i;
  double j;
  double k;
} Quaternion;

std::vector<double> elementGradient(std::vector<int> el, std::vector<double> pts, std::vector<double> phi);
std::vector<double> elementGradient2(std::vector<int> el, std::vector<double> thesePts, std::vector<double> phi);
std::vector<double> pointsInElem(std::vector<int> el, std::vector<double> pts);
std::vector<double> invJacobianTranspose(std::vector<double> thesePts);
int idx(int ix, int jx);

double quaternionDot(const Quaternion a, const Quaternion b);
Quaternion quaternionMult(const Quaternion a, const Quaternion b);
Quaternion quaternionNormalise(const Quaternion q);
void quaternionToRotation(const Quaternion q, double M[3][3]);
Quaternion rotationToQuaternion(double M[3][3]);

double my_dot(const double a[3], const double b[3]);
void my_normalize(double v[3]);
void my_cross(const double a[3], const double b[3], double c[3]);

Quaternion axisSystem(const std::vector<double> apex, std::vector<double> surface_inward);
Quaternion rotateAxis(const Quaternion q, const double alpha, const double beta);
Quaternion slerp(const double a_coeff, const Quaternion a, const double b_coeff, const Quaternion b);
Quaternion biv_interpolate(Quaternion epi_axis, Quaternion lv_axis, Quaternion rv_axis,
                           double lap_epi, double lap_lv, double lap_rv,
                           double alpha_epi, double alpha_endo, double beta_epi, double beta_endo);

void print_vec(const char* name, const double v[3]);
void print_vec(const char* name, const std::vector<double> v);
void print_quaternion(const char* name, const Quaternion q);
void print_matrix(const char* name, double M[3][3]);

int main(int argc, char* argv[]) {
    mitkCommandLineParser parser;

    // Set general information about your command-line app
    parser.setCategory("EASI processing");
    parser.setTitle("Define fibres");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription("Define fibres CemrgApp Command-line App.");

    // How should arguments be prefixed
    parser.setArgumentPrefix("--", "-");

    parser.addArgument(
                "mesh", "msh", mitkCommandLineParser::String,
                "Mesh path.", "Path to mesh (.elem) file.");
    parser.addArgument(
                "alpha_endo", "a_endo", mitkCommandLineParser::String,
                "Alpha_Endo Angle", "angle in degrees");
    parser.addArgument(
                "alpha_epi", "a_epi", mitkCommandLineParser::String,
                "Alpha_Epi Angle", "angle in degrees");
    parser.addArgument(
                "beta_endo", "b_endo", mitkCommandLineParser::String,
                "Beta_Endo Angle", "angle in degrees");
    parser.addArgument(
                "beta_epi", "b_epi", mitkCommandLineParser::String,
                "Beta_Epi Angle", "angle in degrees");
    parser.addArgument(
                "iterations", "iter", mitkCommandLineParser::String,
                "Maximum number of iterations (for debugging)", "Iterations default = 0 (do all elements)");
    parser.addArgument(
                "frequency", "freq", mitkCommandLineParser::String,
                "Prints output at set intervals (for debugging)", "Frequency default = 100000");
    parser.addArgument( // optional
                "verbose", "v", mitkCommandLineParser::Bool,
                "Verbose Output", "Whether to produce verbose output");
    parser.addArgument( // optional
                "debug", "d", mitkCommandLineParser::Bool,
                "Sets debug mode ON", "Debug.");

    // Parse arguments.
    // This method returns a mapping of long argument names to their values.
    std::cout << "Parsing arguments" << '\n';
    auto parsedArgs = parser.parseArguments(argc, argv);

    if (parsedArgs.empty()){
        return EXIT_FAILURE;
    }

    if (parsedArgs["mesh"].Empty()) {
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }

    auto meshfile = us::any_cast<std::string>(parsedArgs["mesh"]);

    // Default values for optional arguments
    std::string alpha_epi_str = "-50";
    std::string alpha_endo_str = "40";
    std::string beta_epi_str = "25";
    std::string beta_endo_str = "-65";
    std::string iterations_str = "0"; // do all elements
    std::string frequency_str = "100000";
    auto debug = false;
    auto verbose = false;

    // Parse, cast and set optional arguments
    std::cout << "Checking for user-defined parameters" << '\n';
    if (parsedArgs.end()!= parsedArgs.find("alpha_epi")){
        alpha_epi_str = us::any_cast<std::string>(parsedArgs["alpha_epi"]);
    }

    if (parsedArgs.end()!= parsedArgs.find("alpha_endo")){
        alpha_endo_str = us::any_cast<std::string>(parsedArgs["alpha_endo"]);
    }

    if (parsedArgs.end()!= parsedArgs.find("beta_epi")){
        beta_epi_str = us::any_cast<std::string>(parsedArgs["beta_epi"]);
    }

    if (parsedArgs.end()!= parsedArgs.find("beta_endo")){
        beta_endo_str = us::any_cast<std::string>(parsedArgs["beta_endo"]);
    }

    if (parsedArgs.end()!= parsedArgs.find("iterations")){
        iterations_str = us::any_cast<std::string>(parsedArgs["iterations"]);
    }

    if (parsedArgs.end()!= parsedArgs.find("frequency")){
        frequency_str = us::any_cast<std::string>(parsedArgs["frequency"]);
    }

    if (parsedArgs.end()!= parsedArgs.find("debug")){
        debug = us::any_cast<bool>(parsedArgs["debug"]);
    }

    if (parsedArgs.end()!= parsedArgs.find("verbose")){
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);
    }


    try{
        MITK_INFO(verbose) << "inputs:";

        MITK_INFO(verbose) << "alpha_epi: " << alpha_epi_str;
        MITK_INFO(verbose) << "alpha_endo: " << alpha_endo_str;
        MITK_INFO(verbose) << "beta_epi: " << beta_epi_str;
        MITK_INFO(verbose) << "beta_endo: " << beta_endo_str;
        MITK_INFO(verbose) << "iterations: " << iterations_str;
        MITK_INFO(verbose) << "frequency: " << frequency_str;

        double alpha_epi = QString::fromStdString(alpha_epi_str).toDouble();
        double alpha_endo = QString::fromStdString(alpha_endo_str).toDouble();
        double beta_epi = QString::fromStdString(beta_epi_str).toDouble();
        double beta_endo = QString::fromStdString(beta_endo_str).toDouble();
        int iterations = QString::fromStdString(iterations_str).toInt();
        int frequency = QString::fromStdString(frequency_str).toInt();
	std::cout << "freq:" << frequency << std::endl;

        QFileInfo fi(QString::fromStdString(meshfile));
        QString path2files = fi.absolutePath() + mitk::IOUtil::GetDirectorySeparator();
        QString basename = fi.baseName();

        QString readInMeshFile = fi.absoluteFilePath();
        QString readInPtsFile = path2files + basename + ".pts";
        QString writeOutFile = path2files + "cemrg_" + basename + "_fibres.lon";

        std::ifstream elemFileRead, ptsFileRead, fapex, flv, frv, fepi;
        std::ofstream outputFileWrite;

        int nElem, nPts;
        MITK_INFO(verbose) << ("[...] Reading in mesh file: " + readInMeshFile).toStdString();
        elemFileRead.open((readInMeshFile).toStdString());
        elemFileRead >> nElem;
        MITK_INFO(verbose) << ("[...] Reading in points file: " + readInPtsFile).toStdString();
        ptsFileRead.open((readInPtsFile).toStdString());
        ptsFileRead >> nPts;

        MITK_INFO(verbose) << ("Loading in " + QString::number(nPts) +" points.").toStdString();
        std::vector<double> pts(3*nPts);
        for(int ix=0; ix<nPts; ix++){
            ptsFileRead >> pts[3*ix + 0];
            ptsFileRead >> pts[3*ix + 1];
            ptsFileRead >> pts[3*ix + 2];
        }
        ptsFileRead.close();

        QString apex_name = path2files + basename + "_psi_ab_potential.dat";
        QString epi_name = path2files + basename + "_phi_epi_potential.dat";
        QString lv_name = path2files + basename + "_phi_lv_potential.dat";
        QString rv_name = path2files + basename + "_phi_rv_potential.dat";

        MITK_INFO(verbose) << "Load in Laplace Solves.";
        std::vector<double> epi_v(nPts, 0.0), apex_v(nPts, 0.0), lv_v(nPts, 0.0), rv_v(nPts, 0.0);

        fapex.open(apex_name.toStdString());
        fepi.open(epi_name.toStdString());
        flv.open(lv_name.toStdString());
        frv.open(rv_name.toStdString());

        for (int ix=0; ix<nPts; ix++){
            fepi >> epi_v[ix];
            fapex >> apex_v[ix];
            flv >> lv_v[ix];
            frv >> rv_v[ix];
        }

        fapex.close();
        fepi.close();
        flv.close();
        frv.close();

        int nIter = (iterations>0) ? iterations : nElem;

        MITK_INFO(verbose) << ("Iterations: " + QString::number(nIter)).toStdString();

        double tol = 1e-7;
        outputFileWrite.open(writeOutFile.toStdString());
        outputFileWrite << "2\n";
        MITK_INFO(verbose) << "Starting rule fibres...";
        for(int qx=0; qx<nIter; qx++){
            std::vector<int> elemIdx(4);
            std::string type;
            int tag;

            elemFileRead >> type;
            elemFileRead >> elemIdx[0];
            elemFileRead >> elemIdx[1];
            elemFileRead >> elemIdx[2];
            elemFileRead >> elemIdx[3];
            elemFileRead >> tag;

            double epi=0.0, lv=0.0, rv=0.0;
            for(int iv=0; iv<4; iv++){
                lv += lv_v[elemIdx[iv]];
                rv += rv_v[elemIdx[iv]];
                epi += epi_v[elemIdx[iv]];
            }
            lv /= 4.0;
            rv /= 4.0;
            epi /= 4.0;

            MITK_INFO(debug) << ("[" + QString::number(qx) + "]").toStdString();

            std::vector<double> epi_gradient, apex_gradient, lv_gradient, rv_gradient;
            std::vector<double> ptsInEl = pointsInElem(elemIdx, pts);

            epi_gradient = elementGradient2(elemIdx, ptsInEl, epi_v);
            apex_gradient = elementGradient2(elemIdx, ptsInEl, apex_v);
            lv_gradient = elementGradient2(elemIdx, ptsInEl, lv_v);
            rv_gradient = elementGradient2(elemIdx, ptsInEl, rv_v);
            if(debug){
                print_vec("epi_gradient" , epi_gradient);
                print_vec("apex_gradient", apex_gradient);
                print_vec("lv_gradient"  , lv_gradient);
                print_vec("rv_gradient"  , rv_gradient);
            }
            std::cout << "print gradients" << '\n';

            Quaternion epi_axis, lv_axis, rv_axis;
            if (epi >= tol){
              epi_axis = axisSystem(apex_gradient, epi_gradient);
              if(debug){
                  print_quaternion("epi_axis", epi_axis);
              }
            }

            if (lv >= tol){
              lv_axis = axisSystem(apex_gradient, lv_gradient);
              if(debug){
                  print_quaternion("lv_axis", lv_axis);
              }
            }

            if (rv >= tol){
              rv_axis = axisSystem(apex_gradient, rv_gradient);
              if(debug){
                  print_quaternion("rv_axis", rv_axis);
              }
            }

            Quaternion q;
            q = biv_interpolate(epi_axis, lv_axis, rv_axis, epi, lv, rv, alpha_epi, alpha_endo, beta_epi, beta_endo);

            double R[3][3];
            quaternionToRotation(q, R);


            // Print out the fiber orientation for this element
            if(qx<26){
                std::cout << qx <<":" << std::setprecision(6) << R[0][0] << " "<< R[1][0] << " "<< R[2][0] << std::endl;
            }
            outputFileWrite << R[0][0] << " "<< R[1][0] << " "<< R[2][0] << " "<< R[0][2] << " "<< R[1][2] << " "<< R[2][2] << '\n';

        }

        elemFileRead.close();
        outputFileWrite.close();

        bool success = true;
        MITK_INFO(success) << "Successful program";
        MITK_INFO(verbose) << "Goodbye!";

    }
    catch (const std::exception &e) {
        MITK_ERROR << e.what();
        return EXIT_FAILURE;
    }
    catch(...) {
        MITK_ERROR << "Unexpected error";
        return EXIT_FAILURE;
    }

}

std::vector<double> pointsInElem(std::vector<int> el, std::vector<double> pts){
    std::vector<double> thesePts(4*3,0);
    for(int ix=0; ix<4; ix++){
        for(int jx=0; jx<3; jx++){
            thesePts[4*ix + jx] = pts[3*el[ix] + jx];
        }
    }
    return thesePts;

}
std::vector<double> elementGradient(std::vector<int> el, std::vector<double> thesePts, std::vector<double> phi){
    std::vector<double> gradient(3,0);
    std::vector<double> thesePhi(4,0);
    std::vector<double> dphi0(4*3,0), dphi(4*3,0);

    dphi0[4*0 + 0] = -1; dphi0[4*0 + 1] = -1; dphi0[4*0 + 2] = -1;
    dphi0[4*1 + 0] =  1; dphi0[4*1 + 1] =  0; dphi0[4*1 + 2] =  0;
    dphi0[4*2 + 0] =  0; dphi0[4*2 + 1] =  1; dphi0[4*2 + 2] =  0;
    dphi0[4*3 + 0] =  0; dphi0[4*3 + 1] =  0; dphi0[4*3 + 2] =  1;

    for(int ix=0; ix<4; ix++){
        thesePhi[ix] = phi[el[ix]];
    }
    // std::cout << "these phi / dphi0" << '\n';
    std::vector<double> invJt = invJacobianTranspose(thesePts);

    for(short int jf=0; jf<4; jf++) {
      for(short int ic=0; ic<3; ic++){
        for(short int jc=0; jc<3; jc++){
          dphi[4*jf + ic]= dphi[4*jf + ic] + invJt[idx(ic,jc)]*dphi0[4*jf+jc];
        }
      }
    }
    // std::cout << "dphi" << '\n';
    // Evaluate the gradient
    for(short int iv=0; iv<4; iv++) {
        for(short int jc=0; jc<3; jc++) {
          gradient[jc] = gradient[jc] + dphi[4*iv + jc] * thesePhi[iv];
        }
    }
    // std::cout << "gradient" << '\n';

    double norm=sqrt(gradient[0]*gradient[0]+gradient[1]*gradient[1]+gradient[2]*gradient[2]);
    if(norm>0.0){
      for(short int jc=0; jc<3; jc++) {
        gradient[jc]=gradient[jc]/norm;
      }
    }
    // std::cout << "norm" << '\n';
    return gradient;
}

std::vector<double> elementGradient2(std::vector<int> el, std::vector<double> thesePts, std::vector<double> phi){
    std::vector<double> gradient(3,0);
    std::vector<double> thesePhi(4,0);
    std::vector<std::vector<double> > dphi0, dphi;

    dphi0.resize(4);
    dphi.resize(4);
    for(short int jf=0; jf<4; jf++){
        dphi0[jf].resize(3);
    }
    for(short int jc=0; jc<3; jc++){
      dphi0[0][jc]=-1;
      dphi0[jc+1][jc]=1;
    }

    for(int ix=0; ix<4; ix++){
        thesePhi[ix] = phi[el[ix]];
    }

    std::vector<double> invJt = invJacobianTranspose(thesePts);

    for(short int jf=0; jf<4; jf++) {
        dphi[jf].resize(3,0);
        for(short int ic=0; ic<3; ic++){
            for(short int jc=0; jc<3; jc++){
                (dphi[jf])[ic]= (dphi[jf])[ic] + invJt[idx(ic,jc)]*(dphi0[jf])[jc];
            }
        }
    }
    // std::cout << "dphi done" << '\n';
    // Evaluate the gradient
    for(short int iv=0; iv<4; iv++) {
        for(short int jc=0; jc<3; jc++) {
          gradient[jc] = gradient[jc] + (dphi[iv])[jc] * thesePhi[iv];
        }
    }
    // std::cout << "gradient done" << '\n';

    double norm=sqrt(gradient[0]*gradient[0]+gradient[1]*gradient[1]+gradient[2]*gradient[2]);
    if(norm>0.0){
      for(short int jc=0; jc<3; jc++) {
        gradient[jc]=gradient[jc]/norm;
      }
    }

    return gradient;

}

double quaternionDot(const Quaternion a, const Quaternion b) {
  return a.s*b.s + a.i*b.i +a.j*b.j +a.k*b.k;
}

Quaternion quaternionMult(const Quaternion a, const Quaternion b) {
  Quaternion result;
  result.s = a.s*b.s - a.i*b.i - a.j*b.j - a.k*b.k;
  result.i = a.s*b.i + a.i*b.s + a.j*b.k - a.k*b.j;
  result.j = a.s*b.j - a.i*b.k + a.j*b.s + a.k*b.i;
  result.k = a.s*b.k + a.i*b.j - a.j*b.i + a.k*b.s;
  return result;
}

Quaternion quaternionNormalise(const Quaternion q) {
  double norm;
  norm = sqrt(quaternionDot(q, q));
  Quaternion result;
  result.s = q.s / norm;
  result.i = q.i / norm;
  result.j = q.j / norm;
  result.k = q.k / norm;
  return result;
}

double my_dot(const double a[3], const double b[3]) {
  double sum=0;
  for (int ii=0; ii<3; ii++) {
    sum += a[ii]*b[ii];
  }
  return sum;
}

void my_normalize(double v[3]) {
  double mag = sqrt(my_dot(v,v));
  for (int ii=0; ii<3; ii++) {
    v[ii] /= mag;
  }
}

void my_cross(const double a[3], const double b[3], double c[3]) {
  c[0] = a[1]*b[2] - b[1]*a[2];
  c[1] = b[0]*a[2] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

Quaternion slerp(const double a_coeff, const Quaternion a, const double b_coeff, const Quaternion b) {
  double tol = 1e-7;

  Quaternion R[4];
  //R0 = a
  R[0] = a;
  //R1 = i*a
  R[1].s = -a.i;
  R[1].i =  a.s;
  R[1].j = -a.k;
  R[1].k =  a.j;
  //R2 = j*a
  R[2].s = -a.j;
  R[2].i =  a.k;
  R[2].j =  a.s;
  R[2].k = -a.i;
  //R3 = k*a
  R[3].s = -a.k;
  R[3].i = -a.j;
  R[3].j =  a.i;
  R[3].k =  a.s;

  //find the quaternion from a canidates that is closest to b.
  double max_dot_p=0;
  int max_quat_index=-1;
  for (int ii=0; ii<4; ii++) {
    double dot_p = quaternionDot(R[ii],b);
    if (dot_p < 0) {
      dot_p = -dot_p;
      R[ii].s *= -1;
      R[ii].i *= -1;
      R[ii].j *= -1;
      R[ii].k *= -1;
    }
    if (dot_p > max_dot_p) {
      max_quat_index=ii;
      max_dot_p = dot_p;
    }
  }
  Quaternion max_quat = R[max_quat_index];

  Quaternion result;
  if (max_dot_p > 1-tol) {
    //in this case, theta == sin_theta == 0.  Quaterions are so close, so just pick one.
    result = b;
  } else {
    //compute the interpolation.
    double theta = acos(max_dot_p);
    double sin_theta = sqrt(1-max_dot_p*max_dot_p);

    double adjusted_a_coeff = sin(theta*a_coeff)/sin_theta;
    double adjusted_b_coeff = sin(theta*b_coeff)/sin_theta;
    result.s = max_quat.s*adjusted_a_coeff + b.s*adjusted_b_coeff;
    result.i = max_quat.i*adjusted_a_coeff + b.i*adjusted_b_coeff;
    result.j = max_quat.j*adjusted_a_coeff + b.j*adjusted_b_coeff;
    result.k = max_quat.k*adjusted_a_coeff + b.k*adjusted_b_coeff;
  }
  return result;
  //return quaternionNormalise(result);
}

void quaternionToRotation(const Quaternion q, double M[3][3]) {
  double w=q.s;
  double x=q.i;
  double y=q.j;
  double z=q.k;
  M[0][0] = w*w+x*x-y*y-z*z;
  M[0][1] = 2*x*y+2*w*z;
  M[0][2] = 2*x*z-2*w*y;

  M[1][0] = 2*x*y-2*w*z;
  M[1][1] = w*w-x*x+y*y-z*z;
  M[1][2] = 2*y*z+2*w*x;

  M[2][0] = 2*x*z+2*w*y;
  M[2][1] = 2*y*z-2*w*x;
  M[2][2] = w*w-x*x-y*y+z*z;
}

Quaternion rotationToQuaternion(double M[3][3]) {
  Quaternion q;
  double tol = 1e-12;

  double w2 = (1+M[0][0]+M[1][1]+M[2][2])/4;
  if (w2 > tol) {
    q.s = sqrt(w2);
    q.i = (M[1][2]-M[2][1])/4/q.s;
    q.j = (M[2][0]-M[0][2])/4/q.s;
    q.k = (M[0][1]-M[1][0])/4/q.s;
  } else {
    q.s = 0;
    double x2 = -(M[1][1]+M[2][2])/2;
    if (x2 > tol) {
      q.i = sqrt(x2);
      q.j = M[0][1]/2/q.i;
      q.k = M[0][2]/2/q.i;
    } else {
      q.i = 0;
      double y2 = (1-M[2][2])/2;
      if (y2 > tol) {
        q.j = sqrt(y2);
        q.k = M[1][2]/2/q.j;
      } else {
        q.j = 0;
        q.k = 1;
      }
    }
  }

  return q;
}

Quaternion axisSystem(const std::vector<double> apex, std::vector<double> surface_inward) {
  //Compute the reference frame
  double M[3][3];

  double e0[3];
  for (int ii=0; ii<3; ii++) {
    e0[ii] = apex[ii];
  }
  my_normalize(e0);
  for (int ii=0; ii<3; ii++) {
    M[ii][0] = e0[ii];
  }
  double e2[3];
  double surf_array[3] = {surface_inward[0], surface_inward[1], surface_inward[2]};
  double dotp = my_dot(e0, surf_array);
  for (int ii=0; ii<3; ii++) {
    e2[ii] = surface_inward[ii] - dotp*e0[ii];
  }
  my_normalize(e2);
  for (int ii=0; ii<3; ii++) {
    M[ii][2] = e2[ii];
  }
  double e1[3];
  my_cross(e2,e0,e1);
  my_normalize(e1);
  for (int ii=0; ii<3; ii++) {
    M[ii][1] = e1[ii];
  }

  return rotationToQuaternion(M);
}

Quaternion rotateAxis(const Quaternion q, const double alpha, const double beta) {
  Quaternion a = {cos(-alpha/2),0,            0, sin(-alpha/2)};
  Quaternion b = {cos(-beta/2) ,sin(-beta/2), 0, 0            };
  return quaternionMult(quaternionMult(b,a),q);
}

Quaternion biv_interpolate(Quaternion epi_axis, Quaternion lv_axis, Quaternion rv_axis,
                           double lap_epi, double lap_lv, double lap_rv,
                           double alpha_epi, double alpha_endo, double beta_epi, double beta_endo) {

  Quaternion q;
  double tol = 1e-7;

  if (lap_lv < tol && lap_rv < tol) {
    // Then we are on the epicardium - just compute epicardial fibres
    q = rotateAxis(epi_axis, (alpha_epi + 90.)*M_PI/180., beta_epi*M_PI/180.);

  } else {

    double lap_septum = lap_lv + lap_rv;
    double alpha_septum = ((90. - alpha_endo) * lap_lv/lap_septum +
                         (90. + alpha_endo) * lap_rv/lap_septum) * M_PI/180.;
    double beta_septum  = (- beta_endo * lap_lv/lap_septum
                         + beta_endo * lap_rv/lap_septum) * M_PI/180.;

    Quaternion lv_q;
    if (lap_lv >= tol) {
      lv_q = rotateAxis(lv_axis, alpha_septum, beta_septum);
    }

    Quaternion rv_q;
    if (lap_rv >= tol) {
      rv_q = rotateAxis(rv_axis, -1.0*alpha_septum, beta_septum);
    }

    Quaternion septum_q;
    if (lap_lv >= tol && lap_rv < tol) {
      septum_q = lv_q;
    } else if (lap_lv < tol && lap_rv >= tol) {
      septum_q = rv_q;
    } else {
      septum_q = slerp(lap_lv/lap_septum, lv_q,
                       lap_rv/lap_septum, rv_q);
    }

    if (lap_epi < tol) {
      q = septum_q;
    } else {
      double alpha_wall = ((alpha_epi  + 90.) * lap_epi +
                         (alpha_endo + 90.) * lap_septum) * M_PI/180.;
      double beta_wall  = (beta_epi * lap_epi
                         + beta_endo * lap_septum) * M_PI/180.;
      Quaternion epi_q = rotateAxis(epi_axis, alpha_wall, beta_wall);
      q = slerp(lap_epi, epi_q, lap_septum, septum_q);
    }
  }

  return q;
}

std::vector<double> invJacobianTranspose(std::vector<double> thesePts){
    std::vector<double> Jt(3*3, 0), invJt(3*3, 0);
    std::vector<double> p0(3,0);
    for (int ix = 0; ix < 3; ix++) {
        p0[ix] = thesePts[4*0 + ix];
    }
    for (int ix = 0; ix < 3; ix++) {
        for (int jx = 0; jx < 3; jx++) {
            Jt[3*ix + jx] = thesePts[4*(ix+1) + jx] - p0[jx];
        }
    }

    //invert the matrix
    double det=0;
    det = Jt[idx(0,0)]*(Jt[idx(1,1)]*Jt[idx(2,2)] - Jt[idx(1,2)]*Jt[idx(2,1)]) +
     -1.0*Jt[idx(0,1)]*( Jt[idx(1,0)]*Jt[idx(2,2)] - Jt[idx(1,2)]*Jt[idx(2,0)] ) +
          Jt[idx(0,2)]*( Jt[idx(1,0)]*Jt[idx(2,1)] - Jt[idx(2,0)]*Jt[idx(1,1)] );

    if(sqrt(det*det)>0.0){
      invJt[idx(0,0)] =        Jt[idx(1,1)]*Jt[idx(2,2)] - Jt[idx(1,2)]*Jt[idx(2,1)];
      invJt[idx(0,1)] = -1.0*( Jt[idx(0,1)]*Jt[idx(2,2)] - Jt[idx(2,1)]*Jt[idx(0,2)]);
      invJt[idx(0,2)] =        Jt[idx(0,1)]*Jt[idx(1,2)] - Jt[idx(0,2)]*Jt[idx(1,1)];

      invJt[idx(1,0)] = -1.0*( Jt[idx(1,0)]*Jt[idx(2,2)] - Jt[idx(1,2)]*Jt[idx(2,0)] );
      invJt[idx(1,1)] =        Jt[idx(0,0)]*Jt[idx(2,2)] - Jt[idx(0,2)]*Jt[idx(2,0)];
      invJt[idx(1,2)] = -1.0*( Jt[idx(0,0)]*Jt[idx(1,2)] - Jt[idx(0,2)]*Jt[idx(1,0)]);

      invJt[idx(2,0)] =        Jt[idx(1,0)]*Jt[idx(2,1)] - Jt[idx(2,0)]*Jt[idx(1,1)];
      invJt[idx(2,1)] = -1.0*( Jt[idx(0,0)]*Jt[idx(2,1)] - Jt[idx(2,0)]*Jt[idx(0,1)]);
      invJt[idx(2,2)] =        Jt[idx(1,1)]*Jt[idx(0,0)] - Jt[idx(1,0)]*Jt[idx(0,1)];

      std::vector<double>::iterator it;
      for(it=invJt.begin(); it!=invJt.end(); ++it){
        *it=*it/det;
      }
    }else {
      invJt.clear();
      invJt.resize(9,0);
    }
    return(invJt);
}

int idx(int ix, int jx){
    return (3*ix+jx);
}

void print_vec(const char* name, const double v[3]) {
  printf("%s = [%g %g %g]\n", name, v[0], v[1], v[2]);
}

void print_vec(const char* name, const std::vector<double> v) {
  printf("%s = [%g %g %g]\n", name, v[0], v[1], v[2]);
}

void print_quaternion(const char* name, const Quaternion q) {
  printf("%s = [%g %g %g %g]\n", name, q.s, q.i, q.j, q.k);
}

void print_matrix(const char* name, double M[3][3]) {
  printf("%s = [\n", name);
  printf(" %g %g %g\n", M[0][0], M[0][1], M[0][2]);
  printf(" %g %g %g\n", M[1][0], M[1][1], M[1][2]);
  printf(" %g %g %g\n", M[2][0], M[2][1], M[2][2]);
  printf("]\n");
}

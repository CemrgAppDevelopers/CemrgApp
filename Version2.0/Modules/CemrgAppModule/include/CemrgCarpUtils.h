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
 *
 * Simple Common Utilities
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef CemrgCarpUtils_h
#define CemrgCarpUtils_h

#include <MitkCemrgAppModuleExports.h>
#include <mitkImage.h>
#include <mitkBoundingObject.h>
#include <mitkDataNode.h>
#include <mitkDataStorage.h>
#include <QString>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

typedef struct {
  double s;
  double i;
  double j;
  double k;
} Quaternion;

class MITKCEMRGAPPMODULE_EXPORT CemrgCarpUtils {
public:
    // basics
    void CemrgCarpUtils();

    inline void SetElementFilename(QSring elementPath){elemPath = elementPath; UpdateFileNames()};
    void UpdateFileNames();
    void ReadElementFile(int nIter=0); // 0 = read whole file
    void ReadPointsFile();
    void ReadLaplaceSolvesFile(QString dataFilePath, int whichFile);
    void ReadElementGradientFile(QString gradientPath, int whichFile);

    inline void SetTolerance(double t){tol = t;};

    inline void SetAlphaEndo(double ae){alpha_endo = ae;};
    inline void SetAlphaEpi(double aE){alpha_epi = aE;};
    inline void SetBetaEndo(double be){beta_endo = be;};
    inline void SetBetaEpi(double bE){beta_epi = bE;};

    inline double AlphaEndo(){return alpha_endo;};
    inline double AlphaEpi(){return alpha_epi;};
    inline double BetaEndo(){return beta_endo;};
    inline double BetaEpi(){return beta_epi;};

    inline int idx(int ix, int jx, int sz=3){return (sz*ix+jx);};

    double DataCentre(std::vector<double> fullDatVector, std::vector<int> elemIdx);
    std::vector<double> GradientAtElement(std::vector<double> fullGradient, std::vector<int> index);

    // Fibres calculation functions
    inline double QuaternionDot(const Quaternion a, const Quaternion b){return (a.s*b.s + a.i*b.i +a.j*b.j +a.k*b.k)};

    Quaternion QuaternionMult(const Quaternion a, const Quaternion b);
    Quaternion QuaternionNormalise(const Quaternion q);
    Quaternion RotationToQuaternion(double M[3][3]);
    void QuaternionToRotation(const Quaternion q, double M[3][3]);

    double DotProd(const double a[3], const double b[3]);
    void NormaliseVector(double v[3]);
    void CrossProd(const double a[3], const double b[3], double c[3]);

    Quaternion AxisSystem(const std::vector<double> apex, std::vector<double> surface_inward);
    Quaternion RotateAxis(const Quaternion q, const double alpha, const double beta);
    Quaternion Slerp(const double a_coeff, const Quaternion a, const double b_coeff, const Quaternion b);
    Quaternion BiVInterpolate(Quaternion epi_axis, Quaternion lv_axis, Quaternion rv_axis, double epi, double lv, double rv);

    void DefineFibres();

    //Static Utils
    static void OriginalCoordinates(QString imagePath, QString pointPath, QString outputPath, double scaling=1000);
    static void CalculateCentreOfGravity(QString pointPath, QString elemPath, QString outputPath);
    static void RegionMapping(QString bpPath, QString pointPath, QString elemPath, QString outputPath);
    static void NormaliseFibreFiles(QString fibresPath, QString outputPath);

    static void RectifyFileValues(QString pathToFile, double minVal=0.0, double maxVal=1.0);

    static int GetTotalFromCarpFile(QString pathToFile, bool totalAtTop=true);
    static std::vector<double> ReadScalarField(QString pathToFile);
    static void CarpToVtk(QString elemPath, QString ptsPath, QString outputPath, bool saveRegionlabels=true);
    static void AppendScalarFieldToVtk(QString vtkPath, QString fieldName, QString typeData, std::vector<double> field, bool setHeader=true);
    static void AppendVectorFieldToVtk(QString vtkPath, QString fieldName, QString typeData, std::vector<double> field, bool setHeader=true);


private:
    std::vector<double> ReadInFile(QString dataFilePath);
    std::vector<double> ReadInGradient(QString gradientPath);

    QString dir, bname, elemPath, ptsPath;
    double alpha_endo, alpha_epi, beta_endo, beta_epi;
    int nPts, nElem;
    double tol;
    std::vector<int> el;
    std::vector<double> pts;
    std::vector<double> apexPts, epiPts, lvPts, rvPts;
    std::vector<double> apexGrad, epiGrad, lvGrad, rvGrad;
};

#endif // CemrgCarpUtils_h

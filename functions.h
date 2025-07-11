#ifndef LILAK_VERSION 
#include "LKCompiled.h"
#include "LKMisc.cpp"
#include "LKPainter.cpp"
#include "LKDrawing.cpp"
#include "LKDrawingGroup.cpp"
#endif

const int kUniformBeamProfile = 0;
const int kGauss2DBeamProfile = 1;
int fBeamTypeList[] = {kGauss2DBeamProfile, kUniformBeamProfile};

TGraph* NewGraph(TString name="graph", int mst=20, double msz=0.6, int mcl=kBlack, int lst=-1, int lsz=-1, int lcl=-1);
TGraphErrors* NewGraphErrors(TString name="graph", int mst=20, double msz=0.6, int mcl=kBlack, int lst=-1, int lsz=-1, int lcl=-1);
double CircleIntersectionArea(double x1, double y1, double r1, double x2, double y2, double r2);
double Gaussian2D(double x, double y, double xc, double yc, double sigmaX, double sigmaY);
double IntegrateGaussian2D(double xc, double yc, double sigmaX, double sigmaY, double x0, double y0, double r0, int nPoints=100);
double IntegrateGaussian2DFast(double xc, double yc, double sigmaX, double sigmaY, double x0, double y0, double r0, int nPoints=100);
double IntegrateGaussian2D_SumGamma(double sigma, double cx, double r0);
double IntegrateGaussian2DFaster(double xc, double yc, double sigmaX, double sigmaY, double x0, double y0, double r0);
const char* BeamTypeString(int beam_type);

class AttenuatorSimulationInput {
    public:
        int constant;
        int exponent;
        int numbering;
        double beam_radius;
        double hole_diameter;
        double attenuation;
        AttenuatorSimulationInput(int f, int e, double d, double r)
        : constant(f), exponent(e), hole_diameter(d), beam_radius(r) {
            attenuation = constant * std::pow(10.0, exponent);
            numbering = constant*100 + (-exponent);
        }
};

double CircleIntersectionArea(double x1, double y1, double r1, double x2, double y2, double r2)
{
    double d = std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    if (d >= r1 + r2) {
        return 0.0;
    }
    if (d <= std::abs(r1 - r2)) {
        double smallerRadius = std::min(r1, r2);
        return M_PI * smallerRadius * smallerRadius;
    }
    double r1Squared = r1 * r1;
    double r2Squared = r2 * r2;
    double angle1 = std::acos((d * d + r1Squared - r2Squared) / (2 * d * r1));
    double angle2 = std::acos((d * d + r2Squared - r1Squared) / (2 * d * r2));
    double triangleArea = 0.5 * std::sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2));
    double segmentArea1 = r1Squared * angle1;
    double segmentArea2 = r2Squared * angle2;
    return segmentArea1 + segmentArea2 - triangleArea;
}

double Gaussian2D(double x, double y, double xc, double yc, double sigmaX, double sigmaY) {
    double norm = 1.0 / (2 * TMath::Pi() * sigmaX * sigmaY);
    double dx = (x - xc) / sigmaX;
    double dy = (y - yc) / sigmaY;
    return norm * TMath::Exp(-0.5 * (dx * dx + dy * dy));
}

double IntegrateGaussian2D(double xc, double yc, double sigmaX, double sigmaY, double x0, double y0, double r0, int nPoints)
{
    double integral = 0.0;
    double dTheta = 2 * TMath::Pi() / nPoints; // Step size in angle
    double dr = r0 / nPoints; // Step size in radius

    for (int i = 0; i < nPoints; ++i) {
        double r = i * dr;
        for (int j = 0; j < nPoints; ++j) {
            double theta = j * dTheta;
            double x = x0 + r * TMath::Cos(theta);
            double y = y0 + r * TMath::Sin(theta);

            // Accumulate the Gaussian value at this point
            integral += Gaussian2D(x, y, xc, yc, sigmaX, sigmaY) * r * dr * dTheta;
        }
    }

    return integral;
}

double IntegrateGaussian2DFast(double xc, double yc, double sigmaX, double sigmaY, double x0, double y0, double r0, int nPoints)
{
    double cx = TMath::Sqrt((x0-xc)*(x0-xc) + (y0-yc)*(y0-yc));
    return IntegrateGaussian2D_SumGamma(sigmaX, cx, r0);
}

double IntegrateGaussian2D_SumGamma(double sigma, double cx, double r0)
{
    cout << "cx = " << cx << endl;
    cx = cx/sigma;
    r0 = r0/sigma;
    double valueSum = 0;
    for (int k=0; k<4; ++k)
    {
        double k_factorial = TMath::Factorial(k);
        double value_add = TMath::Power(cx*cx/2,k)/(k_factorial*k_factorial) * TMath::Gamma(k+1,r0*r0/2);
        cout << k << ") " << TMath::Power(cx*cx/2,k) << " / " << (k_factorial*k_factorial) << " * " << TMath::Gamma(k+1,r0*r0/2) << " = " << value_add << endl;
        valueSum += value_add;
    }
    cout << TMath::Exp(-cx*cx/2) << endl;
    cout << TMath::Exp(-cx*cx/2) * valueSum << endl;
    double integral = 1 - TMath::Exp(-cx*cx/2) * valueSum;
    cout << integral << endl;
    return integral;
}

double IntegrateGaussian2DFaster(double xc, double yc, double sigmaX, double sigmaY, double x0, double y0, double r0)
{
    double dx = (x0 - xc) / sigmaX;
    double dy = (y0 - yc) / sigmaY;
    double gaussianValue = (1.0 / (2 * TMath::Pi() * sigmaX * sigmaY)) *
                           TMath::Exp(-0.5 * (dx * dx + dy * dy));
    double circleArea = TMath::Pi() * r0 * r0;
    return circleArea * gaussianValue;
}

TGraph *NewGraph(TString name, int mst, double msz, int mcl, int lst, int lsz, int lcl)
{
    auto graph = new TGraph();
    graph -> SetName(name);
    if (mst<=0) mst = 20;
    if (msz<=0) msz = 1;
    if (mcl<0) mcl = kBlack;
    graph -> SetMarkerStyle(mst);
    graph -> SetMarkerSize(msz);
    graph -> SetMarkerColor(mcl);
    if (lst<0) lst = 1;
    if (lsz<0) lsz = 1;
    if (lcl<0) lcl = mcl;
    graph -> SetLineStyle(lst);
    graph -> SetLineWidth(lsz);
    graph -> SetLineColor(lcl);
    graph -> SetFillStyle(0);
    return graph;
}

TGraphErrors *NewGraphErrors(TString name, int mst, double msz, int mcl, int lst, int lsz, int lcl)
{
    auto graph = new TGraphErrors();
    graph -> SetName(name);
    if (mst<=0) mst = 20;
    if (msz<=0) msz = 1;
    if (mcl<0) mcl = kBlack;
    graph -> SetMarkerStyle(mst);
    graph -> SetMarkerSize(msz);
    graph -> SetMarkerColor(mcl);
    if (lst<0) lst = 1;
    if (lsz<0) lsz = 1;
    if (lcl<0) lcl = mcl;
    graph -> SetLineStyle(lst);
    graph -> SetLineWidth(lsz);
    graph -> SetLineColor(lcl);
    graph -> SetFillStyle(0);
    return graph;
}

const char* BeamTypeString(int beam_type)
{
    if      (beam_type==kUniformBeamProfile) return "uniform";
    else if (beam_type==kGauss2DBeamProfile) return "gauss2D";
    return "no_beam_type";
}

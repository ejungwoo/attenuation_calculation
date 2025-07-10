#ifndef LILAK_VERSION 
#include "LKCompiled.h"
#include "LKMisc.cpp"
#include "LKPainter.cpp"
#include "LKDrawing.cpp"
#include "LKDrawingGroup.cpp"
#endif
TGraph* NewGraph(TString name="graph", int mst=20, double msz=0.6, int mcl=kBlack, int lst=-1, int lsz=-1, int lcl=-1);
TGraphErrors* NewGraphErrors(TString name="graph", int mst=20, double msz=0.6, int mcl=kBlack, int lst=-1, int lsz=-1, int lcl=-1);
double CircleIntersectionArea(double x1, double y1, double r1, double x2, double y2, double r2);
double Gaussian2D(double x, double y, double xc, double yc, double sigmaX, double sigmaY);
double IntegrateGaussian2D(double xc, double yc, double sigmaX, double sigmaY, double x0, double y0, double r0, int nPoints=100);
double IntegrateGaussian2DFast(double xc, double yc, double sigmaX, double sigmaY, double x0, double y0, double r0, int nPoints=100);
double IntegrateGaussian2D_SumGamma(double sigma, double cx, double r0);
double IntegrateGaussian2DFaster(double xc, double yc, double sigmaX, double sigmaY, double x0, double y0, double r0);
const int kUniformBeamProfile = 0;
const int kGauss2DBeamProfile = 1;
int fBeamTypeList[] = {kGauss2DBeamProfile, kUniformBeamProfile};
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

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
void attenuation_calculation()
{
    bool draw_count_histogram_together = true;

    int constant, exponent;
    double hole_diameter, beam_diameter;
    vector<AttenuatorSimulationInput> attenuator_list = {
        //AttenuatorSimulationInput(constant=1, exponent=-3, hole_diameter=0.12, beam_diameter=5),
        //AttenuatorSimulationInput(constant=1, exponent=-3, hole_diameter=0.12, beam_diameter=10),
        //AttenuatorSimulationInput(constant=5, exponent=-3, hole_diameter=0.12, beam_diameter=5),
        //AttenuatorSimulationInput(constant=5, exponent=-3, hole_diameter=0.12, beam_diameter=10),
        AttenuatorSimulationInput(constant=1, exponent=-5, hole_diameter=0.03, beam_diameter=10),
    };

    /////////////////////////////////////////////////////////////////////

    double attenuator_dx = 40; // width of active attenuator aread (mm)
    double attenuator_half = 0.5*attenuator_dx;
    double attenuator_active_dx = 36; // width of active attenuator area (mm)
    double attenuator_active_dy = 36; // height of active attenuator area (mm)
    double attenuator_radius = 0.5*attenuator_active_dx; // radius of active attenuator area (mm)
    double efficiency_range = 0.30;
    double beam_radius_to_sigma = 1./3;
    int sample_ndivision_x = 200;
    int sample_ndivision_y = 200;
    int rmdr = 1; // ??

    /////////////////////////////////////////////////////////////////////

    auto top = new LKDrawingGroup("attenuator");

    for (auto at : attenuator_list)
    {
        auto drawings = top -> CreateGroup();
        drawings -> SetCanvasSize(500,500,true);

        double beam_sigma = at.beam_radius*beam_radius_to_sigma;
        int beam_size_number = int(at.beam_radius);

        cout << Form("Attenuation = %d*10^{%d},  hole_diameter = %.2f,  beam_radius = %.2f,  resolution range = %.2f",at.constant,at.exponent,at.hole_diameter,at.beam_radius,efficiency_range) << endl;

        auto x0     = 0.5*(attenuator_dx - attenuator_active_dx);
        auto y0     = 0.5*(attenuator_dx - attenuator_active_dy);
        auto area_h = 0.5*at.hole_diameter*0.5*at.hole_diameter*TMath::Pi();
        auto area_a = attenuator_active_dx*attenuator_active_dy;
        auto x_dist = sqrt(2*area_h/sqrt(3)/at.attenuation);
        auto y_dist = sqrt(3)*x_dist/2;
        auto ny     = int(attenuator_active_dy/y_dist)+1;
        auto nx1    = int(attenuator_active_dx/x_dist)+1;
        auto nx     = nx1;

        vector<TVector3> points;
        auto draw_attenuator_sketch = drawings -> CreateDrawing(Form("draw_%d_%d_1_%d",at.numbering,beam_size_number,int(draw_count_histogram_together)));

        TString hname = Form("hist_%dEm%d_h%d_s%d",at.constant,at.exponent,int(100*at.hole_diameter),beam_size_number);
        auto hist = new TH2D(hname,";x (mm);y (mm)",100,0,attenuator_dx,100,0,attenuator_dx);
        hist -> SetStats(0);
        draw_attenuator_sketch -> Add(hist);

        auto graph_holes = NewGraph("graph_holes", 24, (at.exponent<3?0.3:0.4), kBlack);
        auto graph_active_boundary = NewGraph("graph_active_boundary");
        graph_active_boundary -> SetLineStyle(2);
        for (auto i=0; i<=100; ++i)
            graph_active_boundary -> SetPoint(graph_active_boundary->GetN(), attenuator_radius*cos(i*TMath::Pi()/50)+attenuator_half, attenuator_radius*sin(i*TMath::Pi()/50)+attenuator_half);
        draw_attenuator_sketch -> Add(graph_active_boundary,"samel");

        for (auto yy : {y0,attenuator_half,attenuator_dx-y0}) {
            auto line = new TLine(0,yy,attenuator_dx,yy);
            line -> SetLineStyle(2);
            draw_attenuator_sketch -> Add(line,"samel");
        }
        for (auto xx : {x0,attenuator_half,attenuator_dx-x0}) {
            auto line = new TLine(xx,0,xx,attenuator_dx);
            line -> SetLineStyle(2);
            draw_attenuator_sketch -> Add(line,"samel");
        }

        double xLow  = x0;
        double xHigh = x0 + (nx-1)*x_dist;
        double yLow  = y0;
        double yHigh = y0 + (ny-1)*y_dist;
        double dxi = (attenuator_half-0.5*(xLow+xHigh));
        double dyi = (attenuator_half-0.5*(yLow+yHigh));
        //dyi = dyi-0.5*y_dist;
        double xi = x0+dxi;
        double yi = y0+dyi;
        int countY = 0;
        for (auto iy=0; iy<ny; ++iy)
        {
            nx = nx1;
            //if (iy%2==1) nx = nx2;
            double yh = yi + iy*y_dist;
            bool firstX = true;
            TString noteX;
            double xh1, xh2;
            int countX = 0;
            for (auto ix=0; ix<nx; ++ix)
            {
                double xh = xi + ix*x_dist;
                if (iy%2==rmdr) xh = xh + 0.5*x_dist;
                if (xh>x0+attenuator_active_dx) continue;
                if ((attenuator_half-xh)*(attenuator_half-xh)+(attenuator_half-yh)*(attenuator_half-yh)>attenuator_radius*attenuator_radius) continue;
                if (firstX) {
                    firstX = false;
                    noteX = Form(" %d) x=%.4f, y=%.4f",countY+1,xh,yh);
                    xh1 = xh;
                    xh2 = yh;
                }
                graph_holes -> SetPoint(graph_holes->GetN(), xh, yh);
                points.push_back(TVector3(xh,yh,0));
                countX++;
            }
            if (countX>0) countY++;
        }

        draw_attenuator_sketch -> Add(graph_holes,"samep");

        TString title = Form("A=%dx10^{%d}, #phi_{hole}=%.2f, n_{hole}=%d, r_{beam}=%d",at.constant,at.exponent,at.hole_diameter,int(points.size()),beam_size_number);
        hist -> SetTitle(title);

        auto draw_efficiency = drawings -> CreateDrawing(Form("draw_%d_%d_1d_%d",at.numbering,beam_size_number,int(draw_count_histogram_together)),draw_count_histogram_together);
        draw_efficiency -> SetCanvasMargin(0.12,0.12,0.12,0.12);
        draw_efficiency -> SetPaveDx(0.5);
        draw_efficiency -> SetPaveLineDy(0.06);
        draw_efficiency -> SetStatsFillStyle(0);
        draw_efficiency -> SetLegendBelowStats();
        auto lgTG = new TLegend();
        lgTG -> AddEntry((TObject*)0,Form("%dx10^{%d}, #phi=%.2f, r=%d",at.constant,at.exponent,at.hole_diameter,beam_size_number),"");
        lgTG -> SetFillStyle(0);

        vector<TH1D*> hist_efficiency_list;
        for (auto beam_type : fBeamTypeList)
        {
            auto drawAttnRatio = drawings -> CreateDrawing(Form("draw_%d_%d_2d_%s_%d",at.numbering,beam_size_number,BeamTypeString(beam_type),int(draw_count_histogram_together)),draw_count_histogram_together);
            drawAttnRatio -> SetCanvasMargin(0.115,0.14,0.125,0.13);

            auto graphCircle0 = NewGraph("graphCircle0",20,0.6,kBlue);
            auto graphCircleR = NewGraph("graphCircleR"); graphCircleR -> SetLineColor(kBlue);
            auto graphCircle3 = NewGraph("graphCircle3"); graphCircle3 -> SetLineColor(kRed);
            auto graphCircle2 = NewGraph("graphCircle2"); graphCircle2 -> SetLineColor(kRed);
            auto graphCircle1 = NewGraph("graphCircle1"); graphCircle1 -> SetLineColor(kRed);
            auto graphSimCenter = NewGraph("graphSimCenter");
            auto graphSimCenterRange = NewGraph("graphSimCenterRange");
            //double xs1 = attenuator_half - 0.5*x_dist;
            //double xs2 = attenuator_half + 0.5*x_dist;
            //double ys1 = attenuator_half - 0.5*y_dist;
            //double ys2 = attenuator_half + 0.5*y_dist;
            double xs1 = attenuator_half - 1.0*x_dist;
            double xs2 = attenuator_half + 1.0*x_dist;
            double ys1 = attenuator_half - 1.0*y_dist;
            double ys2 = attenuator_half + 1.0*y_dist;
            graphSimCenterRange -> SetPoint(0,xs1,ys1);
            graphSimCenterRange -> SetPoint(1,xs1,ys2);
            graphSimCenterRange -> SetPoint(2,xs2,ys2);
            graphSimCenterRange -> SetPoint(3,xs2,ys1);
            graphSimCenterRange -> SetPoint(4,xs1,ys1);
            int count = 0;
            TString hname = Form("hist_%dEm%d_h%d_s%d",at.constant,at.exponent,int(100*at.hole_diameter),beam_size_number);
            TString title = Form("A=%dx10^{%d}, #phi=%.2f, n_{hole}=%d, r_{beam}=%d",at.constant,at.exponent,at.hole_diameter,int(points.size()),beam_size_number);
            TH1D *hist_efficiency;
            if (beam_type==kUniformBeamProfile) hist_efficiency = new TH1D(hname+"_u_bt"+beam_type+"_1d",title+Form(";Attenuation / %dx10^{%d}",at.constant,at.exponent),100,1-efficiency_range,1+efficiency_range);
            if (beam_type==kGauss2DBeamProfile) hist_efficiency = new TH1D(hname+"_g_bt"+beam_type+"_1d",title+Form(";Attenuation / %dx10^{%d}",at.constant,at.exponent),100,1-efficiency_range,1+efficiency_range);
            hist_efficiency_list.push_back(hist_efficiency);
            hist_efficiency -> SetStats(0);
            if (beam_type==kGauss2DBeamProfile) hist_efficiency -> SetLineColor(kRed);
            auto histSampleCount2 = new TH2D(hname+"_bt"+beam_type+"_2d",title+";beam-center-x (mm); beam-center-y (mm)",sample_ndivision_x,xs1,xs2,sample_ndivision_y,ys1,ys2);
            histSampleCount2 -> SetContour(200);
            histSampleCount2 -> SetStats(0);
            histSampleCount2 -> GetZaxis() -> SetMaxDigits(6);
            histSampleCount2 -> GetZaxis() -> SetDecimals(1);
            for (auto iy=0; iy<sample_ndivision_y; ++iy)
            {
                double y1 = ys1 + iy*(ys2-ys1)/sample_ndivision_y;
                for (auto ix=0; ix<sample_ndivision_x; ++ix)
                {
                    double x1 = xs1 + ix*(xs2-xs1)/sample_ndivision_x;

                    double total = 0;
                    for (auto point : points) {
                        double x2 = point.x();
                        double y2 = point.y();
                        if ( ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)) > at.beam_radius*at.beam_radius)
                            continue;
                        double value = 0;
                        if (beam_type==kUniformBeamProfile) value = CircleIntersectionArea(x1, y1, at.beam_radius, x2, y2, 0.5*at.hole_diameter);
                        if (beam_type==kGauss2DBeamProfile) value = IntegrateGaussian2DFaster(x1, y1, beam_sigma, beam_sigma, x2, y2, 0.5*at.hole_diameter);
                        total += value;
                    }
                    if (beam_type==kUniformBeamProfile) total = total / (TMath::Pi()*at.beam_radius*at.beam_radius);
                    if (beam_type==kGauss2DBeamProfile) total = total / 1;
                    double count_ratio = total/at.attenuation;
                    histSampleCount2 -> SetBinContent(ix+1,iy+1,count_ratio);
                    hist_efficiency -> Fill(count_ratio);
                    //cout << setw(3) << count << setw(10) << x1 << setw(10) << y1 << setw(15) << total << endl;

                    if ((ix==0 && iy==0))
                    {
                        graphCircle0 -> SetPoint(graphCircle0->GetN(),x1,y1);
                        for (auto i=0; i<=100; ++i) {
                            graphCircleR -> SetPoint(graphCircleR->GetN(), at.beam_radius*3*beam_radius_to_sigma*cos(i*TMath::Pi()/50)+x1, at.beam_radius*3*beam_radius_to_sigma*sin(i*TMath::Pi()/50)+y1);
                        }
                    }
                    else if ((ix==sample_ndivision_x-1 && iy==sample_ndivision_y-1))
                    {
                        for (auto i=0; i<=100; ++i) {
                            graphCircle3 -> SetPoint(graphCircle3->GetN(), at.beam_radius*3*beam_radius_to_sigma*cos(i*TMath::Pi()/50)+x1, at.beam_radius*3*beam_radius_to_sigma*sin(i*TMath::Pi()/50)+y1);
                            graphCircle2 -> SetPoint(graphCircle2->GetN(), at.beam_radius*2*beam_radius_to_sigma*cos(i*TMath::Pi()/50)+x1, at.beam_radius*2*beam_radius_to_sigma*sin(i*TMath::Pi()/50)+y1);
                            graphCircle1 -> SetPoint(graphCircle1->GetN(), at.beam_radius*1*beam_radius_to_sigma*cos(i*TMath::Pi()/50)+x1, at.beam_radius*1*beam_radius_to_sigma*sin(i*TMath::Pi()/50)+y1);
                        }
                    }
                    graphSimCenter -> SetPoint(graphSimCenter->GetN(),x1,y1);
                    count++;
                }
            }
            graphSimCenterRange -> SetLineColor(kCyan+1);
            draw_attenuator_sketch -> Add(graphSimCenterRange,"samel");
            graphSimCenter -> SetMarkerStyle(20);
            graphSimCenter -> SetMarkerSize(0.6);
            graphSimCenter -> SetMarkerColor(kRed);
            if (beam_type==kUniformBeamProfile) draw_attenuator_sketch -> Add(graphCircle0,"samel");
            if (beam_type==kUniformBeamProfile) draw_attenuator_sketch -> Add(graphCircleR,"samel");
            if (beam_type==kGauss2DBeamProfile) draw_attenuator_sketch -> Add(graphCircle3,"samel");
            if (beam_type==kGauss2DBeamProfile) draw_attenuator_sketch -> Add(graphCircle2,"samel");
            if (beam_type==kGauss2DBeamProfile) draw_attenuator_sketch -> Add(graphCircle1,"samel");

            drawAttnRatio -> Add(histSampleCount2,"colz");

            auto pv = new TPaveText();
            pv -> SetFillColor(0);
            pv -> SetFillStyle(0);
            pv -> SetBorderSize(0);
            pv -> SetTextAlign(33);
            if (beam_type==kUniformBeamProfile) pv -> AddText("Uniform distribution");
            if (beam_type==kGauss2DBeamProfile) pv -> AddText("2d-Gaus distribution");
            drawAttnRatio -> Add(pv);
            drawAttnRatio -> SetPaveCorner(0);

            if (beam_type==kUniformBeamProfile) lgTG -> AddEntry(hist_efficiency,Form("Uniform (#mu=%.3f,#sigma=%.4f)",hist_efficiency->GetMean(),hist_efficiency->GetStdDev()),"l");
            if (beam_type==kGauss2DBeamProfile) lgTG -> AddEntry(hist_efficiency,Form("2d-Gaus (#mu=%.3f,#sigma=%.4f)",hist_efficiency->GetMean(),hist_efficiency->GetStdDev()),"l");

            auto lg1d = new TLegend();
            lg1d -> AddEntry((TObject*)0,Form("=%dx10^{%d}, #phi=%.2f, r=%d",at.constant,at.exponent,at.hole_diameter,beam_size_number),"");
            lg1d -> SetFillStyle(0);
            if (beam_type==kUniformBeamProfile) lg1d -> AddEntry(hist_efficiency,Form("Uniform (#mu=%.3f,#sigma=%.4f)",hist_efficiency->GetMean(),hist_efficiency->GetStdDev()),"l");
            if (beam_type==kGauss2DBeamProfile) lg1d -> AddEntry(hist_efficiency,Form("2d-Gaus (#mu=%.3f,#sigma=%.4f)",hist_efficiency->GetMean(),hist_efficiency->GetStdDev()),"l");

            auto draw1d = drawings -> CreateDrawing(Form("draw_%d_%d_1d_%s_%d",at.numbering,beam_size_number,BeamTypeString(beam_type),int(draw_count_histogram_together)),!draw_count_histogram_together);
            draw1d -> SetCanvasMargin(0.12,0.12,0.12,0.12);
            draw1d -> SetPaveDx(1);
            draw1d -> SetPaveLineDy(0.05);
            draw1d -> SetStatsFillStyle(0);
            draw1d -> SetLegendBelowStats();
            //cout << hist_efficiency << endl;
            //if (draw1d -> GetEntries()==0)
            draw1d -> Add(hist_efficiency);
            //else
            //    draw1d -> Add(hist_efficiency,"same");
            draw1d -> Add(lg1d);
        }

        if (draw_count_histogram_together) {
            double max = 0;
            for (auto hist_efficiency : hist_efficiency_list) {
                auto max1 = hist_efficiency -> GetMaximum();
                cout << max1 << endl;
                if (max1>max) max = max1;
            }
            for (auto hist_efficiency : hist_efficiency_list) {
                hist_efficiency -> SetMaximum(1.05*max);
                if (draw_efficiency->GetEntries()==0)
                    draw_efficiency -> Add(hist_efficiency);
                else
                    draw_efficiency -> Add(hist_efficiency,"same");
            }
        }
        draw_efficiency -> Add(lgTG);
    }

    top -> Print();
    top -> Draw();
}

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

//double IntegrateGaussian2D(double xc, double yc, double sigmaX, double sigmaY, double x0, double y0, double r0, int nPoints)
//{
//    double integral = 0.0;
//    double dTheta = 2 * TMath::Pi() / nPoints; // Step size in angle
//    double dr = r0 / nPoints; // Step size in radius
//
//    for (int i = 0; i < nPoints; ++i) {
//        double r = i * dr;
//        for (int j = 0; j < nPoints; ++j) {
//            double theta = j * dTheta;
//            double x = x0 + r * TMath::Cos(theta);
//            double y = y0 + r * TMath::Sin(theta);
//
//            // Accumulate the Gaussian value at this point
//            integral += Gaussian2D(x, y, xc, yc, sigmaX, sigmaY) * r * dr * dTheta;
//        }
//    }
//
//    return integral;
//}

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

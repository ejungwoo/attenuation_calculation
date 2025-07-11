#include "functions.h"

void attenuation_calculation()
{
    int constant, exponent;
    double hole_diameter, beam_diameter;
    vector<AttenuatorSimulationInput> attenuator_list = {
        //AttenuatorSimulationInput(constant=1, exponent=-3, hole_diameter=0.12, beam_diameter=5),
        AttenuatorSimulationInput(constant=1, exponent=-3, hole_diameter=0.12, beam_diameter=10),
        //AttenuatorSimulationInput(constant=5, exponent=-3, hole_diameter=0.12, beam_diameter=5),
        AttenuatorSimulationInput(constant=5, exponent=-3, hole_diameter=0.12, beam_diameter=10),
        AttenuatorSimulationInput(constant=1, exponent=-5, hole_diameter=0.03, beam_diameter=10),
    };

    /////////////////////////////////////////////////////////////////////

    bool draw_count_histogram_together = true;
    bool use_larger_beam_center_range = true;

    double attenuator_dx = 40; // width of active attenuator aread (mm)
    double attenuator_active_dx = 36; // width of active attenuator area (mm)
    double attenuator_active_dy = 36; // height of active attenuator area (mm)
    double efficiency_range = 0.30;
    double beam_radius_to_sigma = 1./3;
    int sample_ndivision_x = 200;
    int sample_ndivision_y = 200;
    int grid_type = 0; // 0,1,2

    double attenuator_half = 0.5*attenuator_dx;
    double attenuator_radius = 0.5*attenuator_active_dx; // radius of active attenuator area (mm)

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
        auto hist_attenuator = new TH2D(hname,";x (mm);y (mm)",100,0,attenuator_dx,100,0,attenuator_dx);
        hist_attenuator -> SetStats(0);
        draw_attenuator_sketch -> Add(hist_attenuator);

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
                if (iy%2==grid_type) xh = xh + 0.5*x_dist;
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
        hist_attenuator -> SetTitle(title);

        auto draw_efficiency_all = drawings -> CreateDrawing(Form("draw_%d_%d_1d_%d",at.numbering,beam_size_number,int(draw_count_histogram_together)),draw_count_histogram_together);
        draw_efficiency_all -> SetCanvasMargin(0.12,0.12,0.12,0.12);
        draw_efficiency_all -> SetPaveDx(0.5);
        draw_efficiency_all -> SetPaveLineDy(0.06);
        draw_efficiency_all -> SetStatsFillStyle(0);
        draw_efficiency_all -> SetLegendBelowStats();
        auto lg_efficiency_all = new TLegend();
        lg_efficiency_all -> AddEntry((TObject*)0,Form("%dx10^{%d}, #phi=%.2f, r=%d",at.constant,at.exponent,at.hole_diameter,beam_size_number),"");
        lg_efficiency_all -> SetFillStyle(0);

        vector<TH1D*> hist_efficiency_list;
        for (auto beam_type : fBeamTypeList)
        {
            auto draw_efficiency_2d = drawings -> CreateDrawing(Form("draw_%d_%d_2d_%s_%d",at.numbering,beam_size_number,BeamTypeString(beam_type),int(draw_count_histogram_together)),draw_count_histogram_together);
            draw_efficiency_2d -> SetCanvasMargin(0.115,0.14,0.125,0.13);

            auto graph_circle_0 = NewGraph("graph_circle_0",20,0.6,kBlue);
            auto graph_circle_r = NewGraph("graph_circle_r"); graph_circle_r -> SetLineColor(kBlue);
            auto graph_circle_3 = NewGraph("graph_circle_3"); graph_circle_3 -> SetLineColor(kRed);
            auto graph_circle_2 = NewGraph("graph_circle_2"); graph_circle_2 -> SetLineColor(kRed);
            auto graph_circle_1 = NewGraph("graph_circle_1"); graph_circle_1 -> SetLineColor(kRed);
            auto graph_sim_center = NewGraph("graph_sim_center");
            auto graph_sim_center_range = NewGraph("graph_sim_center_range");
            double xs1 = attenuator_half - 0.5*x_dist;
            double xs2 = attenuator_half + 0.5*x_dist;
            double ys1 = attenuator_half - 0.5*y_dist;
            double ys2 = attenuator_half + 0.5*y_dist;
            if (use_larger_beam_center_range) {
                xs1 = attenuator_half - 1.0*x_dist;
                xs2 = attenuator_half + 1.0*x_dist;
                ys1 = attenuator_half - 1.0*y_dist;
                ys2 = attenuator_half + 1.0*y_dist;
            }
            graph_sim_center_range -> SetPoint(0,xs1,ys1);
            graph_sim_center_range -> SetPoint(1,xs1,ys2);
            graph_sim_center_range -> SetPoint(2,xs2,ys2);
            graph_sim_center_range -> SetPoint(3,xs2,ys1);
            graph_sim_center_range -> SetPoint(4,xs1,ys1);
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
                        graph_circle_0 -> SetPoint(graph_circle_0->GetN(),x1,y1);
                        for (auto i=0; i<=100; ++i) {
                            graph_circle_r -> SetPoint(graph_circle_r->GetN(), at.beam_radius*3*beam_radius_to_sigma*cos(i*TMath::Pi()/50)+x1, at.beam_radius*3*beam_radius_to_sigma*sin(i*TMath::Pi()/50)+y1);
                        }
                    }
                    else if ((ix==sample_ndivision_x-1 && iy==sample_ndivision_y-1))
                    {
                        for (auto i=0; i<=100; ++i) {
                            graph_circle_3 -> SetPoint(graph_circle_3->GetN(), at.beam_radius*3*beam_radius_to_sigma*cos(i*TMath::Pi()/50)+x1, at.beam_radius*3*beam_radius_to_sigma*sin(i*TMath::Pi()/50)+y1);
                            graph_circle_2 -> SetPoint(graph_circle_2->GetN(), at.beam_radius*2*beam_radius_to_sigma*cos(i*TMath::Pi()/50)+x1, at.beam_radius*2*beam_radius_to_sigma*sin(i*TMath::Pi()/50)+y1);
                            graph_circle_1 -> SetPoint(graph_circle_1->GetN(), at.beam_radius*1*beam_radius_to_sigma*cos(i*TMath::Pi()/50)+x1, at.beam_radius*1*beam_radius_to_sigma*sin(i*TMath::Pi()/50)+y1);
                        }
                    }
                    graph_sim_center -> SetPoint(graph_sim_center->GetN(),x1,y1);
                    count++;
                }
            }
            graph_sim_center_range -> SetLineColor(kCyan+1);
            draw_attenuator_sketch -> Add(graph_sim_center_range,"samel");
            graph_sim_center -> SetMarkerStyle(20);
            graph_sim_center -> SetMarkerSize(0.6);
            graph_sim_center -> SetMarkerColor(kRed);
            if (beam_type==kUniformBeamProfile) draw_attenuator_sketch -> Add(graph_circle_0,"samel");
            if (beam_type==kUniformBeamProfile) draw_attenuator_sketch -> Add(graph_circle_r,"samel");
            if (beam_type==kGauss2DBeamProfile) draw_attenuator_sketch -> Add(graph_circle_3,"samel");
            if (beam_type==kGauss2DBeamProfile) draw_attenuator_sketch -> Add(graph_circle_2,"samel");
            if (beam_type==kGauss2DBeamProfile) draw_attenuator_sketch -> Add(graph_circle_1,"samel");

            draw_efficiency_2d -> Add(histSampleCount2,"colz");

            auto pv = new TPaveText();
            pv -> SetFillColor(0);
            pv -> SetFillStyle(0);
            pv -> SetBorderSize(0);
            pv -> SetTextAlign(33);
            if (beam_type==kUniformBeamProfile) pv -> AddText("Uniform distribution");
            if (beam_type==kGauss2DBeamProfile) pv -> AddText("2d-Gaus distribution");
            draw_efficiency_2d -> Add(pv);
            draw_efficiency_2d -> SetPaveCorner(0);

            if (beam_type==kUniformBeamProfile) lg_efficiency_all -> AddEntry(hist_efficiency,Form("Uniform (#mu=%.3f,#sigma=%.4f)",hist_efficiency->GetMean(),hist_efficiency->GetStdDev()),"l");
            if (beam_type==kGauss2DBeamProfile) lg_efficiency_all -> AddEntry(hist_efficiency,Form("2d-Gaus (#mu=%.3f,#sigma=%.4f)",hist_efficiency->GetMean(),hist_efficiency->GetStdDev()),"l");

            auto lg_efficiency_1d = new TLegend();
            lg_efficiency_1d -> AddEntry((TObject*)0,Form("=%dx10^{%d}, #phi=%.2f, r=%d",at.constant,at.exponent,at.hole_diameter,beam_size_number),"");
            lg_efficiency_1d -> SetFillStyle(0);
            if (beam_type==kUniformBeamProfile) lg_efficiency_1d -> AddEntry(hist_efficiency,Form("Uniform (#mu=%.3f,#sigma=%.4f)",hist_efficiency->GetMean(),hist_efficiency->GetStdDev()),"l");
            if (beam_type==kGauss2DBeamProfile) lg_efficiency_1d -> AddEntry(hist_efficiency,Form("2d-Gaus (#mu=%.3f,#sigma=%.4f)",hist_efficiency->GetMean(),hist_efficiency->GetStdDev()),"l");

            auto draw_efficiency_1d = drawings -> CreateDrawing(Form("draw_%d_%d_1d_%s_%d",at.numbering,beam_size_number,BeamTypeString(beam_type),int(draw_count_histogram_together)),!draw_count_histogram_together);
            draw_efficiency_1d -> SetCanvasMargin(0.12,0.12,0.12,0.12);
            draw_efficiency_1d -> SetPaveDx(1);
            draw_efficiency_1d -> SetPaveLineDy(0.05);
            draw_efficiency_1d -> SetStatsFillStyle(0);
            draw_efficiency_1d -> SetLegendBelowStats();
            draw_efficiency_1d -> Add(hist_efficiency);
            draw_efficiency_1d -> Add(lg_efficiency_1d);
        }

        if (draw_count_histogram_together) {
            double max = 0;
            for (auto hist_efficiency : hist_efficiency_list) {
                auto max1 = hist_efficiency -> GetMaximum();
                if (max1>max) max = max1;
            }
            for (auto hist_efficiency : hist_efficiency_list) {
                hist_efficiency -> SetMaximum(1.05*max);
                if (draw_efficiency_all->GetEntries()==0)
                    draw_efficiency_all -> Add(hist_efficiency);
                else
                    draw_efficiency_all -> Add(hist_efficiency,"same");
            }
            draw_efficiency_all -> Add(lg_efficiency_all);
        }
    }

    top -> Draw();
}

#include <iostream>
#include <string>
#include <stdio.h>
//#include <time.h>
#include <sys/stat.h>
#include <TStyle.h>

const std::string currentDateTime() { // Extracts current date and time

    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d___%H-%M-%S", &tstruct);
    return buf; 
    cout << currentDateTime() << endl;
    }

string format_plots(int event_max){

    std::string print_location = "/Users/swathi/Desktop/Final_GP_S/Results/Plots_("+currentDateTime()+")/";
    int result = mkdir(print_location.c_str(), 0777);

    TFile *f = TFile::Open("Simulation_Output.root");
    TTree *t = (TTree*)f->Get("w1w2");

    float line_width = 2;
    gStyle->SetOptStat("emr");
    gStyle->SetStatW(0.15);                
    gStyle->SetStatH(0.1);    
    gStyle->SetStatY(0.88);  
    gStyle->SetStatX(0.88);
    gStyle->SetStatFontSize(0.03);
                           
    //gROOT->SetStyle("Modern");
    gROOT->SetBatch(TRUE);
    //gROOT->ForceStyle();
        
    Int_t color1; 
    Int_t color2;
    color1 = TColor::GetColor("#DAE8FC");
    color2 = TColor::GetColor("#6C8EBF");
    color1 = color2;
     
/**************************** W1 ****************************/

    TCanvas *canvas = new TCanvas("canvas", "", 0, 0, 450, 700);
    
    canvas->Divide(2,3);

    /*------------------------W1PX------------------------*/
    canvas->cd(1);
    TH1D * plot_w1px = new TH1D("plot_w1px","", 80, -400, 400);
    //TCanvas * canvas_w1px = new TCanvas("canvas_w1px", "", 0, 0, 600, 400);

    gPad->SetBottomMargin(0.1);
    plot_w1px->SetFillStyle(3002);
    plot_w1px->SetFillColor(color2);
    plot_w1px->GetXaxis()->SetTitle(" W1 X-Momentum (GeV/c) ");
    plot_w1px->GetYaxis()->SetTitle(" Events ");
    plot_w1px->GetYaxis()->SetRangeUser(0,event_max); // Range
    plot_w1px->Draw("");
    plot_w1px->SetLineWidth(line_width);
    plot_w1px->SetLineColor(color1);

    t->Project("plot_w1px", "w1px_candidate");
    string file_name_1 = "plot_w1px.pdf";
    //canvas_w1px->Print((print_location + file_name_1).c_str());
      
    /*------------------------W1PY------------------------*/
    canvas->cd(3);
    TH1D * plot_w1py = new TH1D("plot_w1py", "", 80, -400, 400);
    //TCanvas * canvas_w1py = new TCanvas("canvas_w1py", "", 0, 0, 600, 400); 

    gPad->SetBottomMargin(0.1);
    plot_w1py->SetFillStyle(3002);
    plot_w1py->SetFillColor(color2);
    plot_w1py->GetXaxis()->SetTitle(" W1 Y-Momentum (GeV/c) ");
    plot_w1py->GetYaxis()->SetTitle(" Events ");
    plot_w1py->GetYaxis()->SetRangeUser(0,event_max); // Range
    plot_w1py->Draw("");
    plot_w1py->SetLineWidth(line_width);
    plot_w1py->SetLineColor(color1);

    t->Project("plot_w1py", "w1py_candidate");
    string file_name_2 = "plot_w1py.pdf";
    //canvas_w1py->Print((print_location + file_name_2).c_str());
   
    /*------------------------W1PZ------------------------*/
    canvas->cd(5);
    TH1D * plot_w1pz = new TH1D("plot_w1pz", "", 80, -400, 400);
    //TCanvas * canvas_w1pz = new TCanvas("canvas_w1pz", "", 0, 0, 600, 400);
    
    gPad->SetBottomMargin(0.1);
    plot_w1pz->SetFillStyle(3002);
    plot_w1pz->SetFillColor(color2);
    plot_w1pz->GetXaxis()->SetTitle(" W1 Z-Momentum (GeV/c) ");
    plot_w1pz->GetYaxis()->SetTitle(" Events ");
    plot_w1pz->GetYaxis()->SetRangeUser(0,event_max); // Range
    plot_w1pz->Draw("");
    plot_w1pz->SetLineWidth(line_width);
    plot_w1pz->SetLineColor(color1);

    t->Project("plot_w1pz", "w1pz_candidate");
    string file_name_3 = "plot_w1pz.pdf";
    //canvas_w1pz->Print((print_location + file_name_3).c_str());

/**************************** W2 ****************************/

    /*------------------------W2PX------------------------*/
    canvas->cd(2);
    TH1D * plot_w2px = new TH1D("plot_w2px", "", 80, -400, 400);
    //TCanvas * canvas_w2px = new TCanvas("canvas_w2px", "", 0, 0, 600, 400);

    gPad->SetBottomMargin(0.1);
    plot_w2px->SetFillStyle(3002);
    plot_w2px->SetFillColor(color2);
    plot_w2px->GetXaxis()->SetTitle(" W2 X-Momentum (GeV/c) ");
    plot_w2px->GetYaxis()->SetTitle(" Events ");
    plot_w2px->GetYaxis()->SetRangeUser(0,event_max); // Range
    plot_w2px->Draw("");
    plot_w2px->SetLineWidth(line_width);
    plot_w2px->SetLineColor(color1);

    t->Project("plot_w2px", "w2px_candidate");
    string file_name_5 = "plot_w2px.pdf";
    //canvas_w2px->Print((print_location + file_name_5).c_str());
      
    /*------------------------W2PY------------------------*/
    canvas->cd(4);
    TH1D * plot_w2py = new TH1D("plot_w2py", "", 80, -400, 400);
    //TCanvas * canvas_w2py = new TCanvas("canvas_w2py", "", 0, 0, 600, 400);

    gPad->SetBottomMargin(0.1);
    plot_w2py->SetFillStyle(3002);
    plot_w2py->SetFillColor(color2); 
    plot_w2py->GetXaxis()->SetTitle(" W2 Y-Momentum (GeV/c) ");
    plot_w2py->GetYaxis()->SetTitle(" Events ");
    plot_w2py->GetYaxis()->SetRangeUser(0,event_max); // Range
    plot_w2py->Draw("");
    plot_w2py->SetLineWidth(line_width);
    plot_w2py->SetLineColor(color1);

    t->Project("plot_w2py", "w2py_candidate");
    string file_name_6 = "plot_w2py.pdf";
    //canvas_w2py->Print((print_location + file_name_6).c_str());
   
    /*------------------------W2PZ------------------------*/
    canvas->cd(6);
    TH1D * plot_w2pz = new TH1D("plot_w2pz", "", 80, -400, 400);
    //TCanvas * canvas_w2pz = new TCanvas("canvas_w1pz", "", 0, 0, 600, 400);

    gPad->SetBottomMargin(0.1);
    plot_w2pz->SetFillStyle(3002);
    plot_w2pz->SetFillColor(color2);
    plot_w2pz->GetXaxis()->SetTitle(" W2 Z-Momentum (GeV/c) ");
    plot_w2pz->GetYaxis()->SetTitle(" Events ");
    plot_w2pz->GetYaxis()->SetRangeUser(0,event_max); // Range
    plot_w2pz->Draw("");
    plot_w2pz->SetLineWidth(line_width);
    plot_w2pz->SetLineColor(color1);

    t->Project("plot_w2pz", "w2pz_candidate");
    string file_name_7 = "plot_w2pz.pdf";
    //canvas_w2pz->Print((print_location + file_name_7).c_str());

    canvas->Print((print_location + "combined_plots.pdf").c_str());

/**************************** DAW ****************************/
    
    TCanvas *canvasd = new TCanvas("canvas", "", 0, 0, 500, 200);
    
    canvasd->Divide(2,1);

    /*------------------------DAW 1------------------------*/
    canvasd->cd(1);
    TH1D * plot_daw1 = new TH1D("plot_daw1", "", 80, 0, 1);
    //TCanvas * canvas_daw1 = new TCanvas("canvas_daw1", "", 0, 0, 600, 400);

    plot_daw1->SetTitle("Angle reconstruction difference for W1");
    plot_daw1->GetXaxis()->SetTitle(" Difference (Rad) ");
    plot_daw1->GetYaxis()->SetTitle(" Events ");
    plot_daw1->GetYaxis()->SetRangeUser(0,2*event_max); // Range
    plot_daw1->SetFillColor(color2);
    plot_daw1->Draw("");
    plot_daw1->SetLineWidth(line_width);
    plot_daw1->SetLineColor(color1);

    t->Project("plot_daw1", "daw1");
    string file_name_4 = "plot_daw1.pdf";
    //canvas_daw1->Print((print_location + file_name_4).c_str());

    mean1.push_back(plot_daw1->GetMean());
    sd1.push_back(plot_daw1->GetStdDev());

    /*------------------------DAW 2------------------------*/
    canvasd->cd(2);
    TH1D * plot_daw2 = new TH1D("plot_daw2", "", 80, 0, 1);
    //TCanvas * canvas_daw2 = new TCanvas("canvas_daw2", "", 0, 0, 600, 400);

    plot_daw2->SetTitle("Angle reconstruction difference for W2");
    plot_daw2->GetXaxis()->SetTitle(" Difference (Rad) ");
    plot_daw2->GetYaxis()->SetTitle(" Events ");
    plot_daw2->GetYaxis()->SetRangeUser(0,2*event_max); // Range
    plot_daw2->SetFillColor(color2);
    plot_daw2->Draw("");
    plot_daw2->SetLineWidth(line_width);
    plot_daw2->SetLineColor(color1);

    t->Project("plot_daw2", "daw2");
    string file_name_8 = "plot_daw2.pdf";

    canvasd->Print((print_location + "combined_daw_plots.pdf").c_str());

    mean2.push_back(plot_daw2->GetMean());
    sd2.push_back(plot_daw2->GetStdDev());

    //FINISHING//
    cout << "***************************************************************************" << endl; 
	cout << "Plotting generation successful." << endl;
    cout << "***************************************************************************" << endl; 

    return print_location;
    }

int format_plots_averages(string print_location, const std::vector<float>& x, const std::vector<float>& s, const std::vector<float>& time_out){

    float line_width = 3;              
    gROOT->SetStyle("Modern");
    gROOT->SetBatch(TRUE);
    gROOT->ForceStyle();

    Int_t color1; 
    Int_t color2;
    color1 = TColor::GetColor("#DAE8FC");
    color2 = TColor::GetColor("#6C8EBF");

    //TCanvas *canvas = new TCanvas("canvas", "", 0, 0, 900, 300);

    //canvas->Divide(3,1);

    string parameter_name = "Energy limit";
    
    //canvas->cd(1);
    TCanvas *canvas1 = new TCanvas("canvas1", "", 0, 0, 180, 150);
    TGraph* graph_daw1 = new TGraph(x.size(), &x[0], &mean1[0]);
    TGraph* graph_daw2 = new TGraph(x.size(), &x[0], &mean2[0]);
    gPad->SetLeftMargin(0.2);
    graph_daw1->SetTitle(""); // Average reconstruction accuracy
    graph_daw1->GetXaxis()->SetTitle(parameter_name.c_str());
    graph_daw1->GetYaxis()->SetTitle("Mean of angle differences (Rad)");
    graph_daw1->SetLineWidth(line_width);
    graph_daw1->SetLineColor(color2);
    graph_daw1->SetMarkerColor(color2);
    graph_daw1->SetMarkerStyle(21);
    graph_daw1->SetMarkerSize(0.3);
    graph_daw1->GetYaxis()->SetRangeUser(0,0.3); // Range
    graph_daw1->Draw("");

    gPad->SetLeftMargin(0.2);
    graph_daw2->SetLineWidth(line_width);
    graph_daw2->SetLineColor(kBlack);
    graph_daw2->SetMarkerColor(kBlack);
    graph_daw2->SetMarkerStyle(21);
    graph_daw2->SetMarkerSize(0.3);
    graph_daw2->GetYaxis()->SetRangeUser(0,0.5); // Range
    graph_daw2->Draw("PL");

    auto legend_p = new TLegend(0.8, 0.8, 0.9, 0.9);
    legend_p->AddEntry(graph_daw1, "W1", "l");
    legend_p->AddEntry(graph_daw2, "W2", "l");
    legend_p->SetEntrySeparation(0.1);
    legend_p->SetMargin(0.5);
    legend_p->Draw();
    canvas1->Print((print_location + "graph_mean.pdf").c_str());

    /*********************************************************/

    //canvas->cd(1);
    TCanvas *canvassd = new TCanvas("canvassd", "", 0, 0, 180, 150);
    TGraph* graph_daw1sd = new TGraph(x.size(), &x[0], &sd1[0]);
    TGraph* graph_daw2sd = new TGraph(x.size(), &x[0], &sd2[0]);
    gPad->SetLeftMargin(0.2);
    graph_daw1sd->SetTitle(""); // Average reconstruction accuracy
    graph_daw1sd->GetXaxis()->SetTitle(parameter_name.c_str());
    graph_daw1sd->GetYaxis()->SetTitle("Standard deviation of angle differences (Rad)");
    graph_daw1sd->SetLineWidth(line_width);
    graph_daw1sd->SetLineColor(color2);
    graph_daw1sd->SetMarkerColor(color2);
    graph_daw1sd->SetMarkerStyle(21);
    graph_daw1sd->SetMarkerSize(0.3);
    graph_daw1sd->GetYaxis()->SetRangeUser(0,0.5); // Range
    graph_daw1sd->Draw("");

    gPad->SetLeftMargin(0.2);
    graph_daw2sd->SetLineWidth(line_width);
    graph_daw2sd->SetLineColor(kBlack);
    graph_daw2sd->SetMarkerColor(kBlack);
    graph_daw2sd->SetMarkerStyle(21);
    graph_daw2sd->SetMarkerSize(0.3);
    graph_daw2sd->GetYaxis()->SetRangeUser(0,0.3); // Range
    graph_daw2sd->Draw("PL");

    auto legend_sd = new TLegend(0.8, 0.8, 0.9, 0.9);
    legend_sd->AddEntry(graph_daw1sd, "W1", "l");
    legend_sd->AddEntry(graph_daw2sd, "W2", "l");
    legend_sd->SetEntrySeparation(0.1);
    legend_sd->SetMargin(0.5);
    legend_sd->Draw();
    canvassd->Print((print_location + "graph_sd.pdf").c_str());

    /*********************************************************/

    //canvas->cd(2);
    TCanvas *canvas2 = new TCanvas("canvas2", "", 0, 0, 150, 150);
    TGraph* graph_daw3 = new TGraph(x.size(), &x[0], &time_out[0]);
    gPad->SetLeftMargin(0.2);
    graph_daw3->SetTitle(""); // Effects on efficiency
    graph_daw3->GetXaxis()->SetTitle(parameter_name.c_str());
    graph_daw3->GetYaxis()->SetTitle("Run time (ms)");
    graph_daw3->SetLineWidth(line_width);
    graph_daw3->SetLineColor(color2);
    graph_daw3->SetMarkerStyle(21);
    graph_daw3->SetMarkerSize(0.3);
    graph_daw3->Draw(""); 
    canvas2->Print((print_location + "graph_runtime.pdf").c_str());

    //canvas->cd(3);
    TCanvas *canvas3 = new TCanvas("canvas3", "", 0, 0, 150, 150);
    TGraph* graph_daw4 = new TGraph(x.size(), &x[0], &s[0]);
    gPad->SetLeftMargin(0.2);
    graph_daw4->SetTitle(""); // Effects on event success
    graph_daw4->GetXaxis()->SetTitle(parameter_name.c_str());
    graph_daw4->GetYaxis()->SetTitle("Success rate (%)");
    graph_daw4->SetLineWidth(line_width);
    graph_daw4->SetLineColor(color2);
    graph_daw4->SetMarkerStyle(21);
    graph_daw4->SetMarkerSize(0.3);
    graph_daw4->Draw(""); 
    canvas3->Print((print_location + "graph_success.pdf").c_str());

    //canvas->Print((print_location + "graph_plots.pdf").c_str());

    return 0;  
    }

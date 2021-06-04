using LorentzVector = ROOT::Math::PxPyPzEVector;

typedef struct {LorentzVector x, p;} StatePoint;

const double c = 300 * 1000 * 1000; // m/s
const double Tesla = c; // eV/m
const double second = c; // m

namespace unit {
  const double n = 0.001 * 0.001 * 0.001;
  const double micro = 0.001 * 0.001;
  const double m = 0.001;
  const double k = 1000;
  const double M = 1000 * 1000;
  const double G = 1000 * 1000 * 1000;
}

const double mass_e = 511 * unit::k; // eV
const double charge_e = -1; // e = 1
const double momentum = 1 * unit::M; // eV
const double b = 1 * unit::m * Tesla; // eV/m
const double radius_expected = (momentum / unit::G) / 0.3 / (b / Tesla); // m


const double width = 2, height = 2; // m
const double padding = 2; // m
const double dx = 0.01, dy = 0.01; // m

const double tau_final = 8 * unit::n * second; // m
const double dtau = 0.01 * unit::n * second; // m

const int fit_initial_n = 200;
const int fit_final_n = 500;

TGraph2D* magnetic_field;

void build_default_magnetic_field(){
  magnetic_field = new TGraph2D();
  const double magnetic_x = width / 2, magnetic_y = height / 2;
  const double padding_x = padding + magnetic_x, padding_y = padding + magnetic_y;
  int i = 0;
  for(double x = - padding_x; x <= padding_x; x += dx){
    for(double y = - padding_y; y <= padding_y; y += dy){
      if(- magnetic_x <= x && x <= magnetic_x && - magnetic_y <= y && y <= magnetic_y){
        magnetic_field->SetPoint(i, x, y, b);
      }else{
        magnetic_field->SetPoint(i, x, y, 0);
      }
      i++;
    }
  }
}

LorentzVector lorentz_force(LorentzVector x, LorentzVector u){
  auto dpdtau = charge_e * magnetic_field->Interpolate(x.X(), x.Y()) * LorentzVector(u.Y(), -u.X(), 0, 0);
  return dpdtau;
}

StatePoint step_next(StatePoint state_point){
  auto x = state_point.x;
  auto p = state_point.p;

  auto dx1 = p / mass_e * dtau;
  auto dp1 = lorentz_force(x, p / mass_e) * dtau;

  auto dx2 = (p + dp1 / 2) / mass_e * dtau;
  auto dp2 = lorentz_force(x + dx1 / 2, (p + dp1 / 2) / mass_e) * dtau;

  auto dx3 = (p + dp2 / 2) / mass_e * dtau;
  auto dp3 = lorentz_force(x + dx2 / 2, (p + dp2 / 2) / mass_e) * dtau;

  auto dx4 = (p + dp3) / mass_e * dtau;
  auto dp4 = lorentz_force(x + dx3, (p + dp3) / mass_e) * dtau;

  auto x_next = x + (dx1 + 2 * dx2 + 2 * dx3 + dx4) / 6;
  auto p_next = p + (dp1 + 2 * dp2 + 2 * dp3 + dp4) / 6;

  StatePoint state_point_next = (StatePoint){x_next, p_next};

  return state_point_next;
}

void runge_kutta(){
  build_default_magnetic_field();
  std::cout << "magnetic field: " << b / c << " Tesla" << std::endl;
  const double energy = TMath::Sqrt(mass_e * mass_e + momentum * momentum);

  const auto initial_coordinates = LorentzVector(-2, 0, 0, 0); // metre, second
  const auto initial_momentum = LorentzVector(momentum, 0, 0, energy); // eV/c, eV
  std::cout << "beta: " << initial_momentum.Beta() << std::endl;
  std::cout << "gamma: " << initial_momentum.Gamma() << std::endl;
  std::cout << "expected length: " << initial_momentum.Beta() * initial_momentum.Gamma() * tau_final;

  StatePoint state_point = (StatePoint){initial_coordinates, initial_momentum};
  
  auto orbit = new TGraph();
  auto fit_initial = new TGraph(1);
  auto fit_final = new TGraph(1);
  auto graph_for_fit = new TGraph();
  orbit->SetPoint(0, state_point.x.X(), state_point.x.Y());

  double tau = 0;
  for(int tau_step = 1; tau < tau_final; tau_step++, tau += dtau){
    state_point = step_next(state_point);
    orbit->SetPoint(tau_step, state_point.x.X(), state_point.x.Y());
    if(tau_step == fit_initial_n){
      fit_initial->SetPoint(0, state_point.x.X(), state_point.x.Y());
    }
    if(tau_step == fit_final_n){
      fit_final->SetPoint(0, state_point.x.X(), state_point.x.Y());
    }
    if(fit_initial_n <= tau_step && tau_step <= fit_final_n){
      graph_for_fit->SetPoint(tau_step - fit_initial_n, state_point.x.X(), state_point.x.Y());
    }
  }

  magnetic_field->SetTitle("Charged particle in a magnetic field simulated in RK4;x [m];y [m];B [c eV/m]");
  auto magnetic_field_histogram = magnetic_field->GetHistogram();

  auto c = new TCanvas("c", "c", 0, 0, 900, 800);
  c->SetRightMargin(0.15);

  auto fit_graph = new TF1("fit_graph", "-1 * sqrt([0] ** 2 - (x - [1]) ** 2) + [2]", fit_initial->GetPointX(0), fit_final->GetPointX(0));
  
  magnetic_field_histogram->SetOption("COLZ");
  magnetic_field_histogram->Draw();
  orbit->Draw("SAME");
  fit_graph->SetLineColor(kRed);
  fit_graph->SetParameters(2.994, -0.993, 2.995);
  fit_graph->SetParNames("radius","x_{center}","y_{center}");
  graph_for_fit->Fit("fit_graph", "R");
  graph_for_fit->Draw("SAME");
  fit_initial->SetMarkerColor(kBlue);
  fit_initial->SetMarkerStyle(20);
  fit_initial->Draw("SAMEP");
  fit_final->SetMarkerColor(kRed);
  fit_final->SetMarkerStyle(20);
  fit_final->Draw("SAMEP");

  gStyle->SetStatX(0.53);
  gStyle->SetStatY(0.86);
  gStyle->SetStatH(0.16);
  gStyle->SetStatW(0.20);
  gStyle->SetOptFit();

  auto legend = new TLegend(0.5, 0.15, 0.8, 0.3);
  legend->AddEntry(orbit, "electron orbit", "l");
  legend->AddEntry(fit_initial, "fit from", "lp");
  legend->AddEntry(fit_final, "fit to", "lp");
  legend->AddEntry(fit_graph, "fitter", "l");
  legend->Draw();

  auto latex = new TLatex();
  latex->SetTextSize(0.025);
  latex->DrawLatex(-2, -1.50, Form("p_{i} = %.4e GeV/c", momentum / unit::G));
  latex->DrawLatex(-2, -1.75, Form("B = %.4e Tesla", b / Tesla));
  latex->DrawLatex(-2, -2.25, "p[GeV/c] = 0.3B[T]R[m]");
  latex->DrawLatex(-2, -2.50, Form("R_{expected} = %.4e m", radius_expected));

  TArrow *ar2 = new TArrow(-2, 0.2, -1.5, 0.2, 0.02, "|>");
  ar2->SetAngle(40);
  ar2->SetLineWidth(2);
  ar2->Draw();
}
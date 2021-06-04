using LorentzX = ROOT::Math::XYZTVector;
using LorentzP = ROOT::Math::PxPyPzEVector;
using LorentzP_M = ROOT::Math::PxPyPzMVector;

const double c = 300 * 1000 * 1000; // m/s
const double Tesla = c; // eV/m
const double second = c; // m

namespace unit {
  const double n = 0.001 * 0.001 * 0.001;
  const double micro = 0.001 * 0.001;
  const double m = 0.001;
  const double c = 0.01;
  const double k = 1000;
  const double M = 1000 * 1000;
  const double G = 1000 * 1000 * 1000;
}

const double mass_e = 511 * unit::k; // eV
const double charge_e = -1; // e = 1
const double momentum = 2 * unit::M; // eV

const double tau_final = 1 * unit::n * second; // m
const double dtau = 0.001 * unit::n * second; // m

const double edge_margin = 0.005;

class DrainRectangle;

class BeamRK4 {
  public:
    BeamRK4(const LorentzX, const LorentzP_M, const double, TH2D*, const double, const double);
    void step_RK4();
    void plot_orbit_point();
    void set_magnetic_unit(const double magnetic_unit);
    void set_length_unit(const double length_unit);
    void add_drain_rectangle(DrainRectangle*);
    bool is_anihilated();

    LorentzX x;
    LorentzP p;

    const double charge;
    
    TH2D* magnetic_field;
    double magnetic_field_unit = 1.0; //c*eV/m
    double length_unit = 1.0; // m

    vector<DrainRectangle*> drain_rectangles;

    TGraph* orbit;
    int tau_index = 0;
    const double dtau;
    const double tau_final;

  private:
    LorentzX get_dx(const LorentzP p);
    LorentzP get_dp(const LorentzX x, const LorentzP p);
};

class DrainRectangle{
  public:
    DrainRectangle(const double x1, const double x2, const double y1, const double y2)
      : x1(x1), x2(x2), y1(y1), y2(y2) {}
    virtual bool is_collided(BeamRK4* beam_RK4) {
      return (x1 <= beam_RK4->x.X() && beam_RK4->x.X() <= x2 && y1 <= beam_RK4->x.Y() && beam_RK4->x.Y() <= y2);
    }
    void set_length_unit(const double length_unit) {
      this->length_unit = length_unit;
      tbox = new TBox(x1 / length_unit, y1 / length_unit, x2 / length_unit, y2 / length_unit);
      tbox->SetFillStyle(3001);
      tbox->SetFillColor(kRed);
    }

    const double x1, x2, y1, y2;
    double length_unit;
    TBox* tbox;
};

class Detector : public DrainRectangle {
  public:
    Detector(const double x1, const double x2, const double y1, const double y2) : DrainRectangle(x1, x2, y1, y2) {
      hist = new TH1D("hist", "hist", 500, 0, 3 * unit::M);
    }
    bool is_collided(BeamRK4* beam_RK4) override {
      bool ret = (x1 <= beam_RK4->x.X() && beam_RK4->x.X() <= x2 && y1 <= beam_RK4->x.Y() && beam_RK4->x.Y() <= y2);
      if(ret)
        hist->Fill(beam_RK4->p.E());
      return ret;
    }
    TH1D* hist;
};

BeamRK4::BeamRK4(const LorentzX initial_x, const LorentzP_M initial_p, const double particle_charge, TH2D* _magnetic_field, const double dtau, const double tau_final)
: x(initial_x), p(initial_p), charge(particle_charge), magnetic_field(_magnetic_field), dtau(dtau), tau_final(tau_final){
  orbit = new TGraph();
}

void BeamRK4::plot_orbit_point(){
  orbit->SetPoint(tau_index, x.X() / length_unit, x.Y() / length_unit);
  tau_index++;
}

LorentzX BeamRK4::get_dx(const LorentzP p){
  return p / this->p.M() * dtau;
}

LorentzP BeamRK4::get_dp(const LorentzX x, const LorentzP p){
  auto dp = charge * magnetic_field->Interpolate(x.X() / length_unit, x.Y() / length_unit) * magnetic_field_unit * LorentzP(p.Y(), -p.X(), 0, 0) / this->p.M() * dtau;
  return dp;
}

void BeamRK4::set_magnetic_unit(const double magnetic_field_unit){
  this->magnetic_field_unit = magnetic_field_unit;
}

void BeamRK4::set_length_unit(const double length_unit){
  this->length_unit = length_unit;
}

void BeamRK4::add_drain_rectangle(DrainRectangle* drain_rectangle){
  drain_rectangles.emplace_back(drain_rectangle);
}

bool BeamRK4::is_anihilated(){
  if(tau_index * dtau > tau_final){
    // std::cout << "anihilated by time limit" << std::endl;
    return true;
  }
  if(x.X() - edge_margin <= magnetic_field->GetXaxis()->GetXmin() * length_unit || magnetic_field->GetXaxis()->GetXmax() * length_unit <= x.X() + edge_margin || x.Y() - edge_margin <= magnetic_field->GetYaxis()->GetXmin() * length_unit || magnetic_field->GetYaxis()->GetXmax() * length_unit <= x.Y() + edge_margin){
    // std::cout << "anihilated by colliding with outer boundary" << std::endl;
    return true;
  }
  if(std::any_of(drain_rectangles.begin(), drain_rectangles.end(), [this](DrainRectangle* drain_rectangle){return drain_rectangle->is_collided(this);})){
    // std::cout << "anihilated by colliding with rectangle" << std::endl;
    return true;
  }
  return false;
}

void BeamRK4::step_RK4(){
  auto dx1 = get_dx(p);
  auto dp1 = get_dp(x, p);

  auto dx2 = get_dx(p + dp1 / 2);
  auto dp2 = get_dp(x + dx1 / 2, p + dp1 / 2);

  auto dx3 = get_dx(p + dp2 / 2);
  auto dp3 = get_dp(x + dx2 / 2, p + dp2 / 2);

  auto dx4 = get_dx(p + dp3);
  auto dp4 = get_dp(x + dx3, p + dp3);

  x += (dx1 + 2 * dx2 + 2 * dx3 + dx4) / 6;
  p += (dp1 + 2 * dp2 + 2 * dp3 + dp4) / 6;
}

void runge_kutta_spectro(){
  TCanvas* c1 = new TCanvas("c1", "test");

  double cm = unit::c;

  TFile* file = TFile::Open("mfield.root");
  auto MagneticField = file->Get<TH2D>("MagneticField"); //mTesla
  MagneticField->Draw("COLZ");

  DrainRectangle* top_collimator_left = new DrainRectangle(-2 * cm, -1.1 * cm, 3 * cm, 4 * cm);
  top_collimator_left->set_length_unit(cm);
  top_collimator_left->tbox->Draw();
  DrainRectangle* top_collimator_right = new DrainRectangle(-0.9 * cm, 0 * cm, 3 * cm, 4 * cm);
  top_collimator_right->set_length_unit(cm);
  top_collimator_right->tbox->Draw();

  DrainRectangle* side_collimator_top = new DrainRectangle(3 * cm, 4 * cm, -0.9 * cm, 0 * cm);
  side_collimator_top->set_length_unit(cm);
  side_collimator_top->tbox->Draw();
  DrainRectangle* side_collimator_bottom = new DrainRectangle(3 * cm, 4 * cm, -2 * cm, -1.1 * cm);
  side_collimator_bottom->set_length_unit(cm);
  side_collimator_bottom->tbox->Draw();

  Detector* detector = new Detector(4.5 * cm, 5 * cm, -1.5 * cm, -0.5 * cm);
  detector->set_length_unit(cm);
  detector->tbox->SetFillColor(kGreen);
  detector->tbox->Draw();

  TRandom* trandom_momentum = new TRandom();
  TRandom* trandom_angle = new TRandom();

  for(int i = 0; i < 2000000; i++){
    const auto initial_coordinates = LorentzX(-1 * unit::c, 4 * unit::c, 0, 0); // metre, second

    double momentum_amount = trandom_momentum->Gaus(momentum, momentum);
    double angle = trandom_angle->Rndm() * 2 * TMath::Pi(); 
    const auto initial_momentum = LorentzP_M(momentum_amount * TMath::Cos(angle), momentum_amount * TMath::Sin(angle), 0, mass_e); // eV/c, eV

    BeamRK4 beam_RK4 = BeamRK4(initial_coordinates, initial_momentum, charge_e, MagneticField, dtau, tau_final);
    beam_RK4.set_magnetic_unit(unit::m * Tesla);
    beam_RK4.set_length_unit(unit::c);
    beam_RK4.add_drain_rectangle(top_collimator_left);
    beam_RK4.add_drain_rectangle(top_collimator_right);
    beam_RK4.add_drain_rectangle(side_collimator_top);
    beam_RK4.add_drain_rectangle(side_collimator_bottom);
    beam_RK4.add_drain_rectangle(detector);

    beam_RK4.plot_orbit_point();

    while(!beam_RK4.is_anihilated()){
      beam_RK4.step_RK4();
      beam_RK4.plot_orbit_point();
    }

    // std::cout << "anihilated at (" << beam_RK4.x.X() << ", " << beam_RK4.x.Y() << ")" << std::endl; 

    beam_RK4.orbit->Draw("SAME");
  }
  c1->SaveAs("test.png");

  TCanvas* c2 = new TCanvas("c2", "test");
  detector->hist->Draw();
  c2->SaveAs("test_hist.png");
}
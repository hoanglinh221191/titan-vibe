#include "titan_electrostatics.h"
#include "titan_gb_model.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <span>

namespace {

void require(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL titan_gb_model_smoke: " << message << '\n';
    std::exit(1);
  }
}

bool near(double observed, double expected, double tolerance = 1.0e-12) {
  return std::abs(observed - expected) <=
         tolerance * std::max(1.0, std::abs(expected));
}

} // namespace

int main() {
  require(near(titan_gb::screening_factor_for_element('H'), 0.85),
          "hydrogen screen factor");
  require(near(titan_gb::screening_factor_for_element('C'), 0.72),
          "carbon screen factor");
  require(near(titan_gb::screening_factor_for_element('S'), 0.96),
          "sulfur screen factor");

  // Independent numeric oracle for a carbon target (r=1.70 A) descreened by
  // one oxygen source (r=1.50 A, scale=0.85) at 2.50 A.  Values were evaluated
  // from the published/reference HCT analytical pair expression, not generated
  // by the implementation under test.
  const double offset_carbon = 1.70 - titan_gb::kDielectricOffsetAngstrom;
  const double scaled_oxygen =
      (1.50 - titan_gb::kDielectricOffsetAngstrom) * 0.85;
  const double integral = titan_gb::hct_pair_integral_inverse_angstrom(
      offset_carbon, scaled_oxygen, 2.50);
  require(near(integral, 0.01694340732981706), "HCT pair integral oracle");
  require(near(titan_gb::born_radius_angstrom(1.70, integral, 1),
               1.6551506660013076),
          "HCT Born radius oracle");
  require(near(titan_gb::born_radius_angstrom(1.70, integral, 2),
               1.644065677468608),
          "OBC1 Born radius oracle");
  require(near(titan_gb::born_radius_angstrom(1.70, integral, 5),
               1.6518887474844162),
          "OBC2 Born radius oracle");
  require(near(titan_gb::born_radius_angstrom(1.70, 0.0, 5), 1.61),
          "isolated atom returns offset radius");
  require(std::isnan(titan_gb::hct_pair_integral_inverse_angstrom(
              offset_carbon, scaled_oxygen, 0.0)),
          "coincident atoms fail closed");

  const double kappa = titan_electrostatics::amber_gb_kappa_nm_inverse(
      0.15, 78.5, 298.15);
  require(near(kappa, 0.9301946053189923), "Amber GB kappa oracle");
  require(near(std::exp(-kappa), 0.3944769355915672),
          "Amber GB screening oracle");
  require(near(titan_electrostatics::amber_gb_kappa_nm_inverse(
                   0.0, 78.5, 298.15),
               0.0),
          "zero salt has zero kappa");

  const double f_gb =
      titan_gb::effective_gb_distance_nm(0.16, 0.18, 0.22);
  require(near(f_gb, 0.41763813857057064),
          "effective GB pair distance oracle");
  const double pair_energy = titan_gb::polar_pair_energy_kj_mol(
      0.5, -0.25, 0.16, 0.18, 0.22, kappa, 1.0, 78.5);
  require(near(pair_energy, 41.224483479311886),
          "screened GB polar pair oracle");
  require(near(titan_gb::polar_pair_energy_kj_mol(
                   0.5, -0.25, 0.16, 0.18, 0.22, kappa, 1.0, 78.5,
                   false, 0.001),
               pair_energy) &&
              near(titan_gb::polar_pair_energy_kj_mol(
                       0.5, -0.25, 0.16, 0.18, 0.22, kappa, 1.0, 78.5,
                       false),
                   pair_energy),
          "ALPB-off must preserve the established GB pair result");
  const double self_energy = 0.5 * titan_gb::polar_pair_energy_kj_mol(
                                       0.5, 0.5, 0.0, 0.18, 0.18, kappa,
                                       1.0, 78.5);
  require(near(self_energy, -95.443360261148015),
          "screened GB polar self oracle");
  const double intramolecular_endpoint_delta =
      titan_gb::polar_pair_energy_kj_mol(0.5, -0.25, 0.16, 0.20, 0.24,
                                         kappa, 1.0, 78.5) -
      pair_energy;
  require(near(intramolecular_endpoint_delta, -0.73800093441225556),
          "free-to-complex intramolecular pair delta oracle");

  // Independent Sigalov ALPB Eq. 25 oracle.  alpha is dimensionless and A is
  // the whole-system electrostatic size, not an empirical length shift.
  const double alpb_pair = titan_gb::polar_pair_energy_kj_mol(
      0.5, -0.25, 0.16, 0.18, 0.22, kappa, 1.0, 78.5, true, 1.7);
  require(near(alpb_pair, 40.999760876098684),
          "Sigalov ALPB screened pair oracle");
  const double alpb_self =
      0.5 * titan_gb::polar_pair_energy_kj_mol(
                0.5, 0.5, 0.0, 0.18, 0.18, kappa, 1.0, 78.5, true, 1.2);
  require(near(alpb_self, -94.857094530019893),
          "Sigalov ALPB screened self oracle");
  require(std::isnan(titan_gb::polar_pair_energy_kj_mol(
              0.5, -0.25, 0.16, 0.18, 0.22, kappa, 1.0, 78.5, true, 0.0)),
          "ALPB must fail closed without a positive electrostatic size");

  const std::array<double, 1> sphere_x = {11.0};
  const std::array<double, 1> sphere_y = {-7.0};
  const std::array<double, 1> sphere_z = {3.0};
  const std::array<double, 1> sphere_radius = {2.0};
  require(near(titan_gb::electrostatic_size_adet_nm(
                   sphere_x, sphere_y, sphere_z, sphere_radius),
               0.2),
          "ADet electrostatic size of one sphere equals its radius");

  const std::array<double, 3> shape_x = {0.0, 2.0, -1.0};
  const std::array<double, 3> shape_y = {0.0, 1.0, 3.0};
  const std::array<double, 3> shape_z = {0.0, -2.0, 1.0};
  const std::array<double, 3> shape_radius = {1.5, 1.7, 1.2};
  const double shape_size = titan_gb::electrostatic_size_adet_nm(
      shape_x, shape_y, shape_z, shape_radius);
  require(near(shape_size, 0.2895266211543685),
          "ADet asymmetric-shape oracle");
  const std::array<double, 3> translated_x = {4.0, 6.0, 3.0};
  const std::array<double, 3> translated_y = {-5.0, -4.0, -2.0};
  const std::array<double, 3> translated_z = {7.0, 5.0, 8.0};
  require(near(titan_gb::electrostatic_size_adet_nm(
                   translated_x, translated_y, translated_z, shape_radius),
               shape_size),
          "ADet electrostatic size must be translation invariant");

  // Full complex - receptor - ligand ALPB endpoint-cycle oracle.  The three
  // ADet sizes come from distinct endpoint geometries and atomic radii; the
  // expected totals were evaluated independently from Sigalov Eqs. 7, 23, 25.
  const std::array<double, 4> endpoint_x_a = {0.0, 3.5, 8.0, 10.5};
  const std::array<double, 4> endpoint_y_a = {0.0, 0.0, 1.0, -1.0};
  const std::array<double, 4> endpoint_z_a = {0.0, 0.0, 0.0, 0.5};
  const std::array<double, 4> endpoint_radii_a = {1.5, 1.7, 1.2, 1.4};
  const double receptor_size = titan_gb::electrostatic_size_adet_nm(
      std::span<const double>(endpoint_x_a.data(), 2),
      std::span<const double>(endpoint_y_a.data(), 2),
      std::span<const double>(endpoint_z_a.data(), 2),
      std::span<const double>(endpoint_radii_a.data(), 2));
  const double ligand_size = titan_gb::electrostatic_size_adet_nm(
      std::span<const double>(endpoint_x_a.data() + 2, 2),
      std::span<const double>(endpoint_y_a.data() + 2, 2),
      std::span<const double>(endpoint_z_a.data() + 2, 2),
      std::span<const double>(endpoint_radii_a.data() + 2, 2));
  const double complex_size = titan_gb::electrostatic_size_adet_nm(
      endpoint_x_a, endpoint_y_a, endpoint_z_a, endpoint_radii_a);
  require(near(receptor_size, 0.25329511674118116) &&
              near(ligand_size, 0.21958347591280861) &&
              near(complex_size, 0.41558000558572822),
          "three-endpoint ADet oracle");

  const std::array<double, 4> endpoint_q = {0.4, -0.3, 0.2, -0.1};
  const std::array<double, 4> endpoint_x_nm = {0.0, 0.35, 0.8, 1.05};
  const std::array<double, 4> endpoint_y_nm = {0.0, 0.0, 0.1, -0.1};
  const std::array<double, 4> endpoint_z_nm = {0.0, 0.0, 0.0, 0.05};
  const std::array<double, 4> complex_born = {0.21, 0.23, 0.18, 0.20};
  const std::array<double, 2> receptor_born = {0.19, 0.22};
  const std::array<double, 2> ligand_born = {0.17, 0.19};
  const auto endpoint_energy =
      [&](std::span<const double> q, std::span<const double> x,
          std::span<const double> y, std::span<const double> z,
          std::span<const double> born, double size_nm) {
        double energy = 0.0;
        for (std::size_t i = 0; i < q.size(); ++i) {
          energy += 0.5 * titan_gb::polar_pair_energy_kj_mol(
                              q[i], q[i], 0.0, born[i], born[i], kappa,
                              1.0, 78.5, true, size_nm);
          for (std::size_t j = i + 1; j < q.size(); ++j) {
            const double dx = x[i] - x[j];
            const double dy = y[i] - y[j];
            const double dz = z[i] - z[j];
            energy += titan_gb::polar_pair_energy_kj_mol(
                q[i], q[j], dx * dx + dy * dy + dz * dz, born[i], born[j],
                kappa, 1.0, 78.5, true, size_nm);
          }
        }
        return energy;
      };
  const double complex_energy =
      endpoint_energy(endpoint_q, endpoint_x_nm, endpoint_y_nm, endpoint_z_nm,
                      complex_born, complex_size);
  const double receptor_energy = endpoint_energy(
      std::span<const double>(endpoint_q.data(), 2),
      std::span<const double>(endpoint_x_nm.data(), 2),
      std::span<const double>(endpoint_y_nm.data(), 2),
      std::span<const double>(endpoint_z_nm.data(), 2), receptor_born,
      receptor_size);
  const double ligand_energy = endpoint_energy(
      std::span<const double>(endpoint_q.data() + 2, 2),
      std::span<const double>(endpoint_x_nm.data() + 2, 2),
      std::span<const double>(endpoint_y_nm.data() + 2, 2),
      std::span<const double>(endpoint_z_nm.data() + 2, 2), ligand_born,
      ligand_size);
  require(near(complex_energy, -43.72335810533599) &&
              near(receptor_energy, -41.955794724396284) &&
              near(ligand_energy, -11.75222273471625) &&
              near(complex_energy - receptor_energy - ligand_energy,
                   9.9846593537765447),
          "three-endpoint ALPB binding-cycle oracle");

#if defined(__AVX__)
  const std::array<double, 4> avx_q_i = {0.5, -0.4, 0.2, 0.7};
  const std::array<double, 4> avx_q_j = {-0.25, 0.3, 0.2, -0.1};
  const std::array<double, 4> avx_r2 = {0.16, 0.09, 0.0, 0.49};
  const std::array<double, 4> avx_born_i = {0.18, 0.21, 0.17, 0.25};
  const std::array<double, 4> avx_born_j = {0.22, 0.19, 0.17, 0.20};
  std::array<double, 4> avx_q_product{};
  std::array<double, 4> avx_f_gb{};
  std::array<double, 4> avx_screen{};
  for (std::size_t lane = 0; lane < avx_q_i.size(); ++lane) {
    avx_q_product[lane] = avx_q_i[lane] * avx_q_j[lane];
    avx_f_gb[lane] = titan_gb::effective_gb_distance_nm(
        avx_r2[lane], avx_born_i[lane], avx_born_j[lane]);
    avx_screen[lane] = std::exp(-kappa * avx_f_gb[lane]);
  }
  const __m256d avx_energy =
      titan_gb::alpb_polar_pair_energy_from_fgb_avx(
          _mm256_loadu_pd(avx_q_product.data()),
          _mm256_loadu_pd(avx_f_gb.data()),
          _mm256_loadu_pd(avx_screen.data()), 1.0, 78.5, complex_size);
  std::array<double, 4> avx_observed{};
  _mm256_storeu_pd(avx_observed.data(), avx_energy);
  for (std::size_t lane = 0; lane < avx_q_i.size(); ++lane) {
    const double scalar = titan_gb::polar_pair_energy_kj_mol(
        avx_q_i[lane], avx_q_j[lane], avx_r2[lane], avx_born_i[lane],
        avx_born_j[lane], kappa, 1.0, 78.5, true, complex_size);
    require(near(avx_observed[lane], scalar),
            "ALPB AVX/scalar lane parity");
  }
#else
  require(false, "ALPB AVX/scalar parity smoke requires AVX support");
#endif

  std::cout << "PASS titan_gb_model_smoke\n";
  return 0;
}

#include "arguments.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

std::string arguments::get_build_number() const {
  return "titan_fmm_kifmm_parser_smoke";
}
int arguments::help_print() { return 0; }

namespace {

namespace fs = std::filesystem;

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_fmm_kifmm_parser_smoke: " << message << '\n';
    std::exit(1);
  }
}

void parse(arguments &args, std::vector<std::string> values) {
  std::vector<char *> argv;
  argv.reserve(values.size());
  for (std::string &value : values) {
    argv.push_back(value.data());
  }
  require(args.parse_arguments(static_cast<int>(argv.size()), argv.data()) == 0,
          "parser returned success");
}

void check_defaults() {
  arguments args;
  require(args.pbe.tabi_fmm_engine == "cartesian",
          "FMM engine defaults to cartesian");
  require(args.pbe.tabi_kifmm_order == 6, "KIFMM order defaults to 6");
  require(args.pbe.tabi_kifmm_equiv_scale == 1.05,
          "KIFMM equivalent scale default");
  require(args.pbe.tabi_kifmm_check_scale == 1.60,
          "KIFMM check scale default");
  require(args.pbe.tabi_kifmm_solver == "qr", "KIFMM solver default");
  require(args.pbe.tabi_kifmm_profile == "off", "KIFMM profile default");
}

void check_parser_contract() {
  arguments lite;
  parse(lite, {"titan", "-pbe-tabi-fmm-engine", "kifmm-lite",
               "-pbe-tabi-kifmm-order", "8",
               "-pbe-tabi-kifmm-equiv-scale", "1.10",
               "-pbe-tabi-kifmm-check-scale", "1.70",
               "-pbe-tabi-kifmm-solver", "svd",
               "-pbe-tabi-kifmm-profile", "summary"});
  require(lite.pbe.tabi_fmm_engine == "kifmm-lite",
          "parser stores kifmm-lite engine");
  require(lite.pbe.tabi_kifmm_order == 8, "parser stores KIFMM order");
  require(lite.pbe.tabi_kifmm_equiv_scale == 1.10,
          "parser stores equivalent scale");
  require(lite.pbe.tabi_kifmm_check_scale == 1.70,
          "parser stores check scale");
  require(lite.pbe.tabi_kifmm_solver == "svd", "parser stores solver");
  require(lite.pbe.tabi_kifmm_profile == "summary",
          "parser stores profile");

  arguments local;
  parse(local, {"titan", "--pbe-tabi-fmm-engine", "kifmm-local",
                "--pbe-tabi-kifmm-order", "4",
                "--pbe-tabi-kifmm-solver", "normal-tikhonov",
                "--pbe-tabi-kifmm-profile", "full"});
  require(local.pbe.tabi_fmm_engine == "kifmm-local",
          "parser stores kifmm-local engine");
  require(local.pbe.tabi_kifmm_order == 4, "parser stores split order");
  require(local.pbe.tabi_kifmm_solver == "normal-tikhonov",
          "parser stores normal-tikhonov solver");
  require(local.pbe.tabi_kifmm_profile == "full",
          "parser stores full profile");

  arguments experimental;
  parse(experimental,
        {"titan", "-pbe-tabi-fmm-engine", "kifmm-lite-experimental",
         "-pbe-tabi-kifmm-order", "6", "-pbe-tabi-kifmm-equiv-scale",
         "1.05", "-pbe-tabi-kifmm-check-scale", "2.40",
         "-pbe-tabi-kifmm-solver", "qr", "-pbe-tabi-kifmm-profile",
         "summary"});
  require(experimental.pbe.tabi_fmm_engine == "kifmm-lite-experimental",
          "parser stores experimental KIFMM GMRES engine");
  require(experimental.pbe.tabi_kifmm_order == 6,
          "parser stores experimental KIFMM order");
  require(experimental.pbe.tabi_kifmm_check_scale == 2.40,
          "parser stores experimental check scale");

  const fs::path arg_dir = fs::temp_directory_path() / "titan_kifmm_args";
  fs::create_directories(arg_dir);
  const fs::path arg_file = arg_dir / "arg.in";
  {
    std::ofstream out(arg_file);
    out << "-pbe-tabi-fmm-engine=cartesian\n";
    out << "-pbe-tabi-kifmm-order=6\n";
    out << "-pbe-tabi-kifmm-equiv-scale=1.05\n";
    out << "-pbe-tabi-kifmm-check-scale=1.60\n";
    out << "-pbe-tabi-kifmm-solver=qr\n";
    out << "-pbe-tabi-kifmm-profile=off\n";
  }
  arguments from_file;
  parse(from_file, {"titan", "-a", arg_file.string()});
  require(from_file.pbe.tabi_fmm_engine == "cartesian",
          "arg.in parser stores cartesian engine");
  require(from_file.pbe.tabi_kifmm_order == 6,
          "arg.in parser stores KIFMM order");
  fs::remove_all(arg_dir);
}

int run_probe(const std::string &probe) {
  arguments args;
  if (probe == "--probe-invalid-engine") {
    parse(args, {"titan", "-pbe-tabi-fmm-engine", "bogus"});
  } else if (probe == "--probe-invalid-solver") {
    parse(args, {"titan", "-pbe-tabi-kifmm-solver", "bad"});
  } else if (probe == "--probe-invalid-profile") {
    parse(args, {"titan", "-pbe-tabi-kifmm-profile", "debug"});
  }
  return 0;
}

void expect_reject(const fs::path &executable, const std::string &probe,
                   const std::string &expected) {
  const fs::path log = executable.parent_path() / (probe.substr(2) + ".log");
  const std::string command =
      "\"\"" + executable.string() + "\" " + probe + " > \"" +
      log.string() + "\" 2>&1\"";
  require(std::system(command.c_str()) != 0, probe + " exits with failure");
  std::ifstream in(log);
  const std::string text((std::istreambuf_iterator<char>(in)),
                         std::istreambuf_iterator<char>());
  require(text.find(expected) != std::string::npos,
          probe + " emits expected validation message");
}

void check_rejections(const fs::path &executable) {
  expect_reject(executable, "--probe-invalid-engine",
                "-pbe-tabi-fmm-engine");
  expect_reject(executable, "--probe-invalid-solver",
                "-pbe-tabi-kifmm-solver");
  expect_reject(executable, "--probe-invalid-profile",
                "-pbe-tabi-kifmm-profile");
}

} // namespace

int main(int argc, char **argv) {
  if (argc > 1) {
    return run_probe(argv[1]);
  }
  check_defaults();
  check_parser_contract();
  check_rejections(fs::absolute(argv[0]));
  std::cout << "PASS titan_fmm_kifmm_parser_smoke\n";
  return 0;
}

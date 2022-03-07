import * as wasm from "boltzmann";

const S0 = 2891.2;
const RHOC = 1.05375e-5;

document.getElementById("btnsolve").addEventListener("click", solve, false);

// document.getElementById("n").addEventListener("input", solve, false);
// document.getElementById("sigma").addEventListener("input", solve, false);
// document.getElementById("x-start").addEventListener("input", solve, false);
// document.getElementById("x-end").addEventListener("input", solve, false);
// document.getElementById("mass").addEventListener("input", solve, false);

function checkWasmSupported() {
  try {
    if (
      typeof WebAssembly === "object" &&
      typeof WebAssembly.instantiate === "function"
    ) {
      const module = new WebAssembly.Module(
        Uint8Array.of(0x0, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00)
      );
      if (module instanceof WebAssembly.Module)
        return new WebAssembly.Instance(module) instanceof WebAssembly.Instance;
    }
  } catch (e) {}
  return false;
}

function alertNeg(name, value) {
  alert("Invalid " + name + " = " + value + ". Must be > 0.");
}

function getConfig() {
  const config = wasm.Config.new();
  var n = parseInt(document.getElementById("n").value);
  var sigma = parseFloat(document.getElementById("sigma").value);
  var x_start = parseFloat(document.getElementById("x-start").value);
  var x_end = parseFloat(document.getElementById("x-end").value);
  var mass = parseFloat(document.getElementById("mass").value);

  if (n < 0) {
    alertNeg("n", n);
    n = 0;
  }
  if (mass <= 0) {
    alertNeg("mass", mass);
    mass = 100.0;
  }
  if (x_start <= 0) {
    alertNeg("x-start", mass);
    x_start = 1.0;
  }
  if (x_end <= 0) {
    alertNeg("x-end", mass);
    x_end = 500.0;
  }
  if (sigma <= 0) {
    alertNeg("sigma", sigma);
    sigma = 1e-9;
  }
  if (x_end <= x_start) {
    alert(
      "Invalid start and end x = (" +
        x_start +
        ", " +
        x_end +
        ")" +
        ". End must be greater than start."
    );
    x_start = 1;
    x_end = 500;
  }

  config.n = n;
  config.sigma = sigma;
  config.x_start = x_start;
  config.x_end = x_end;
  config.mass = mass;

  return config;
}

function plot(solution, mass) {
  const xs = [];
  const ys = [];
  // console.log("Making plot.");

  if (solution.get_status() != "Success") {
    alert(solution.get_status());
    return;
  }

  const size = solution.size();
  // console.log("size = " + size);
  var wmin = Infinity;
  var wmax = -Infinity;
  for (let i = 0; i < size; ++i) {
    const w = solution.get_w(i);
    const x = Math.exp(solution.get_logx(i));
    const y = Math.exp(w);
    if (w > wmax) {
      wmax = w;
    }
    if (w < wmin) {
      wmin = w;
    }
    xs.push(x);
    ys.push(y);
    // console.log("x, Y = " + x + " " + y);
  }

  const omegah2 = (S0 / RHOC) * mass * ys[ys.length - 1];

  wmin = Math.floor(Math.log10(Math.exp(wmin))) - 1;
  wmax = Math.ceil(Math.log10(Math.exp(wmax))) + 1;
  const num_ws = wmax - wmin;
  const ws = [];
  for (let i = 0; i < num_ws; ++i) {
    ws.push(Math.pow(10, i + wmin));
  }
  console.log(ws[0]);
  console.log(ws[num_ws - 1]);

  const data = [
    {
      type: "scatter",
      x: xs,
      y: ys,
    },
  ];

  const layout = {
    title: {
      text: "$\\Omega h^2 = " + omegah2 + "$",
    },
    xaxis: {
      title: {
        text: "$\\huge{x = m_{\\chi} / T}$",
      },
      type: "log",
      autorange: true,
      linewidth: 2,
      showline: true,
      mirror: "ticks",
    },
    yaxis: {
      title: {
        text: "$\\Large{Y = n_{\\chi}(T)/s(T)}$",
        font: {
          size: 20,
        },
      },
      linewidth: 2,
      mirror: "ticks",
      type: "log",
      // autorange: true,
      tickvals: ws,
      tickformat: "1.0e",
      nticks: 10,
    },
    height: 640,
  };
  const other = { displayModeBar: false };

  Plotly.newPlot("plotlydiv", data, layout, other);
}

function solve() {
  const solution = wasm.solve(getConfig());
  const mass = parseFloat(document.getElementById("mass").value);
  plot(solution, mass);
}

solve();

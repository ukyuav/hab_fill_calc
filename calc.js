var inputs = {
    c_d : null,
    vel : null,
    M_he : null,  // HELIUM MOLAR MASS (kg/mol)
    M_air : null, // AIR MOLAR MASS (kg/mol)
    R : null,// GAS CONSTANT (J* mol^(-1) * K^(-1))
    g : null, //ACCELERATION DUE TO GRAVITY (m*s^(-2))
    V_tank : null, // Internal volume of tank (in^3)
    // Sorry about English units on this one but we need it for psi
    m_balloon_g : null,
    m_payload_g : null,
    p_mb_input : null,
    temp_c_input : null,
}

function getInputValue(){
    inputs.c_d = Number(document.getElementById("c_d_input").value);
    inputs.vel = Number(document.getElementById("vel_input").value);
    inputs.M_he = Number(document.getElementById("M_he_input").value);  
    inputs.M_air = Number(document.getElementById("M_air_input").value); 
    inputs.R = Number(document.getElementById("R_input").value);
    inputs.g = Number(document.getElementById("g_input").value); 
    inputs.V_tank = Number(document.getElementById("V_tank_input").value);
    inputs.m_balloon_g = Number(document.getElementById("m_balloon_g_input").value);
    inputs.m_payload_g = Number(document.getElementById("m_payload_g_input").value);
    inputs.p_mb_input = Number(document.getElementById("p_mb_input").value);
    inputs.temp_c_input = Number(document.getElementById("temp_c_input").value);

    console.log(inputs.c_d);
    console.log(inputs.vel);
    console.log(inputs.M_he);
    console.log(inputs.M_air);
    console.log(inputs.R);
    console.log(inputs.g);
    console.log(inputs.V_tank);
    console.log(inputs.m_balloon_g);
    console.log(inputs.m_payload_g);
    console.log(inputs.p_mb_input);
    console.log(inputs.temp_c_input);
}

function calculate(){
    getInputValue();
    console.log(inputs.c_d);
    m_balloon = inputs.m_balloon_g/1000.0;
    console.log(m_balloon);
    m_payload = inputs.m_payload_g/1000.0;
    console.log(m_payload);
    p_air = inputs.p_mb_input * 100.0;
    console.log(p_air);
    temp = 273.15 + inputs.temp_c_input;
    console.log(temp);

    var rho_he = (p_air * inputs.M_he)/(inputs.R * temp); // Density of Helium (kg/m^3)
    var rho_air = (p_air * inputs.M_air)/(inputs.R * temp); // Density of Air (kg/m^3)
    console.log(rho_he);
    console.log(rho_air);

    var buoyant_accel = ((rho_air - rho_he) * inputs.g) / rho_he // m*s^-2
    var drag_accel = inputs.c_d * 0.5  * rho_air * Math.pow(inputs.vel,2) * Math.PI * Math.pow((3.0/(4.0 * Math.PI))/ rho_he,2.0/3.0) // m*s^-2
    var mass_load = m_balloon + m_payload // kg

    var a = (buoyant_accel/inputs.g) - 1.0; // Unitless
    var b = -drag_accel/inputs.g; // Unitless
    var c =  - mass_load; // kg

    // Expanded solution for ax + bx^(2/3) + c = 0
    // Symbolic solution was obtained with Sympy's utilities.codegen.jscode function
    // The numeric solution is used here reduce response time for use on webpage.
    // The solution, mass_he_needed, is in kilograms.
    var mass_he_needed = -0.111111111111111*(-9.0*Math.pow(c, 2)/Math.pow(a, 2) + 9.0*Math.pow(Math.pow(a, 2)*c + 0.333333333333333*Math.pow(b, 3), 2)/Math.pow(a, 6))/Math.cbrt(Math.sqrt(-Math.pow(-Math.pow(c, 2)/Math.pow(a, 2) + Math.pow(Math.pow(a, 2)*c + 0.333333333333333*Math.pow(b, 3), 2)/Math.pow(a, 6), 3) + Math.pow(0.5*Math.pow(c, 3)/Math.pow(a, 3) - 0.5*Math.pow(c, 2)*(3.0*Math.pow(a, 2)*c + Math.pow(b, 3))/Math.pow(a, 5) + Math.pow(Math.pow(a, 2)*c + 0.333333333333333*Math.pow(b, 3), 3)/Math.pow(a, 9), 2)) + 0.5*Math.pow(c, 3)/Math.pow(a, 3) - 0.5*Math.pow(c, 2)*(3.0*Math.pow(a, 2)*c + Math.pow(b, 3))/Math.pow(a, 5) + Math.pow(Math.pow(a, 2)*c + 0.333333333333333*Math.pow(b, 3), 3)/Math.pow(a, 9)) - 1.0*Math.cbrt(Math.sqrt(-Math.pow(-Math.pow(c, 2)/Math.pow(a, 2) + Math.pow(Math.pow(a, 2)*c + 0.333333333333333*Math.pow(b, 3), 2)/Math.pow(a, 6), 3) + Math.pow(0.5*Math.pow(c, 3)/Math.pow(a, 3) - 0.5*Math.pow(c, 2)*(3.0*Math.pow(a, 2)*c + Math.pow(b, 3))/Math.pow(a, 5) + Math.pow(Math.pow(a, 2)*c + 0.333333333333333*Math.pow(b, 3), 3)/Math.pow(a, 9), 2)) + 0.5*Math.pow(c, 3)/Math.pow(a, 3) - 0.5*Math.pow(c, 2)*(3.0*Math.pow(a, 2)*c + Math.pow(b, 3))/Math.pow(a, 5) + Math.pow(Math.pow(a, 2)*c + 0.333333333333333*Math.pow(b, 3), 3)/Math.pow(a, 9)) - 0.333333333333333*(3.0*Math.pow(a, 2)*c + Math.pow(b, 3))/Math.pow(a, 3)

    var V_m = (mass_he_needed * inputs.R * temp)/(inputs.M_he * p_air) // Volume of gas needed (m^3)

    var V_i = (V_m/0.0283168)*(12**3) // Volume of Gas needed (in^3)

    var P_i = p_air/6895; // Atmospheric pressure in psi (lb*in^-2)
    var Pf = (P_i/inputs.V_tank) * V_i;

    var gross_lift = (buoyant_accel * mass_he_needed)/inputs.g;

    var free_lift = gross_lift - (m_payload + m_balloon);

    console.log(mass_he_needed);
}

/* var c_d = 0.3;
var vel = 5.0;
var M_he = 0.0040026;  // HELIUM MOLAR MASS (kg/mol)
var M_air = 0.02896; // AIR MOLAR MASS (kg/mol)
var R = 8.314; // GAS CONSTANT (J* mol^(-1) * K^(-1))
var g = 9.81; //ACCELERATION DUE TO GRAVITY (m*s^(-2))
var V_tank = 2990 // Internal volume of tank (in^3)
// Sorry about English units on this one but we need it for psi

var m_balloon = 0.6;
var m_payload = 0.1;
var p_air = 101300.0;
var temp = 273.15 + 32;
 */
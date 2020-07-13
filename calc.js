var inputs = {
    c_d : null,
    vel : null,
    M_he : null,  // HELIUM MOLAR MASS (kg/mol)
    M_air : null, // AIR MOLAR MASS (kg/mol)
    R : null,// GAS CONSTANT (J* mol^(-1) * K^(-1))
    g : null, //ACCELERATION DUE TO GRAVITY (m*s^(-2))
    V_tank : null, // Internal volume of tank (in^3)
    // Sorry about English units on this one but we need it for psi
    m_balloon_g : null, // Mass of balloon (g)
    m_payload_g : null, // Mass of payload (g)
    p_mb_input : null, // Surface air pressure (millibars)
    temp_c_input : null, // Surface air temp (degrees Celsius)
}

var outputs = {
    mass : null, // Mass of gas needed (kg)
    volume_l: null, // Volume of Gas needed (l)
    volume_m3: null, // Volume of Gas needed (m^3)
    volume_ft3: null, // Volume of Gas needed (ft^3)
    volume_in3: null, // Volume of Gas needed (in^3)
    psi: null, // Pressure needed to be released from tank (psi)
    gross_lift: null, // Lift created by the gas (kg)
    free_lift: null, // Amount of lift beyond what is needed to support the mass of the balloon and payload (kg)
}
const error_messages = {
    1 : 'All inputs except for temperature should be positive.',
    2 : 'The molar mass of your lifting gas must be less than that of air.'
}

function validate(){
    var errors = [];

    //Check for problems with input

    //Check that all inputs (except temp) are positive
    for(const property in inputs){
        if(property != "temp_c_input" && inputs[property] < 0){
            errors.push(1);
        }
    }

    //Check that we are using a lighter-than-air-gas
    if(inputs.M_he >= inputs.M_air){
        errors.push(2);
    }

    showError(errors);
    return errors.length;
}

function showError(errors){
    var alert_area = document.getElementById("alert-area");
    var alert = document.getElementById("alert-0");

    // Clear previous errors
    while (alert_area.lastChild.id !== 'alert-0') {
        alert_area.removeChild(alert_area.lastChild);
    }

    if(errors.length ===0){
        alert_area.hidden = true;
    }
    else{ // There's a problem with the data
        // Hide the output
        document.getElementById("results-area").hidden = true;
        //  Put up an alert
        alert_area.hidden = false;
        alert.innerHTML = "Please fix missing or incorrect fields.";

        for(error in errors){
            var alert_id = "alert-" + errors[error];
            alert_area.insertAdjacentHTML("beforeend", '<div class="alert alert-danger" role="alert"></div>');
            var new_alert = alert_area.lastElementChild;
            new_alert.setAttribute("id", alert_id);
            new_alert.innerHTML = error_messages[errors[error]];
        }
    }
}

function outputAnswer(mass_he_needed){
    var results_area = document.getElementById("results-area");
    results_area.hidden = false;
    for(const property in outputs){
        document.getElementById(property).innerHTML = outputs[property].toFixed(3);
    }
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
}

function calculate(){    
    getInputValue();

    var invalid = validate();
    if(invalid != 0){
        return false;
    }

    var m_balloon = inputs.m_balloon_g/1000.0; // Mass of balloon in kg
    var m_payload = inputs.m_payload_g/1000.0; // Mass of payload in kg
    var p_air = inputs.p_mb_input * 100.0; // Surface air pressure in pascals
    var temp = 273.15 + inputs.temp_c_input; // Temp in Kelvin
    var rho_he = (p_air * inputs.M_he)/(inputs.R * temp); // Density of Helium (kg/m^3)
    var rho_air = (p_air * inputs.M_air)/(inputs.R * temp); // Density of Air (kg/m^3)

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

    var gross_lift = (buoyant_accel * mass_he_needed)/inputs.g; // Amount of lift created by the gas

    outputs.mass        = mass_he_needed;
    outputs.volume_m3   = V_m;
    outputs.volume_in3  = V_i;
    outputs.volume_ft3  = V_i /(12**3); 
    outputs.volume_l    = V_m * 1000;
    outputs.psi         = (P_i/inputs.V_tank) * V_i;
    outputs.gross_lift  = gross_lift;
    outputs.free_lift   = gross_lift - (m_payload + m_balloon);

    outputAnswer();
    return false;
}
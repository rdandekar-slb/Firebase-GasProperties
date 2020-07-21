function main()
{
    document.getElementById("gas-properties").innerHTML="";
    input = getparams();
    if (input.spgr=="" || input.pres=="" || input.temp=="")
    {
        //alert("Please check data");
        setTimeout(() => {
            document.getElementById("gas-properties").innerHTML="";
        }, 1000);
        document.getElementById("gas-properties").innerHTML="<h3>Please check data</h3>";

        return;
    }
    gasprops=getgasproperties(input);
    strGasProps=JSON.stringify(gasprops);
    document.getElementById("gas-properties").innerHTML=strGasProps;
}
function getparams()
{
    spgr=document.getElementById("spgr").value;
    //alert(spgr);
    pres=document.getElementById("pres").value;
    temp=document.getElementById("temp").value;
    var inputObj = {spgr:spgr, pres:pres, temp:temp};
    return inputObj;
}
function getgasproperties(inp)
{
    h2s = inp.h2s ? inp.h2s : 0.00;
    co2 = inp.co2 ? inp.co2 : 0.00;
    var pres=Number(inp.pres);
    var temp=Number(inp.temp);
    var spgr = parseFloat(inp.spgr);


    r = 10.73;
    a1 = 0.3265;
    a2 = -1.07;
    a3 = -0.5339;
    a4 = 0.01569;
    a5 = -0.05165;
    a6 = 0.5475;
    a7 = -0.7361;
    a8 = 0.1844;
    a9 = 0.1056;
    a10 = 0.6134;
    a11 = 0.721;

    epsilon = (120.0 * (Math.pow((h2s + co2),0.9) - Math.pow((h2s + co2),1.6)));
    epsilon += (15.0 * (Math.pow(h2s,0.5) - Math.pow(h2s,4)));
    ppc = 756.8 - (131.0 * spgr) - (3.6 * spgr * spgr);
    tpc = 169.2 + (349.5 * spgr) - (74.0 * spgr * spgr);

    //console.log(epsilon, ppc, tpc, temp, pres, spgr, h2s, co2);

    tpc -= epsilon;
    ppc = ppc * tpc / (tpc + h2s * (1 - h2s) * epsilon);
    ppr = pres / ppc;
    tpr = (temp + 459.67) / tpc;

    z_temp = 1.0;

    while (true)
    {
        rhopr=0.27*(ppr/(z_temp*tpr));
        c1=a1+(a2/tpr)+(a3/tpr**3)+(a4/tpr**4)+(a5/tpr**5);
        c2=a6+(a7/tpr)+(a8/tpr**2);
        c3=a9*((a7/tpr)+(a8/tpr**2));
        c4=a10*(1+a11*rhopr**2)*((rhopr**2/tpr**3)*Math.exp(-a11*rhopr**2));

        fz=z_temp - (1+ c1*rhopr + c2*rhopr**2 - c3*rhopr**5 + c4);
        dfz= 1 + (c1*rhopr/z_temp) + (2*c2*rhopr**2/z_temp) - (5*c3*rhopr**5/z_temp)+ (2 * (a10*rhopr**2) * ((1 + a11*rhopr**2-(a11*rhopr**2)**2) / (tpr**3 * z_temp)) * Math.exp(-a11 * rhopr ** 2));
        z_new = z_temp - fz/dfz;
        if (Math.abs(z_new - z_temp) < 1e-6) 
        {
            break;

        }
        else
        {
            z_temp = z_new;
        }
        
    }    
    zfact = z_new;
    b_g = 0.02827*zfact*(temp+459.67)/pres;
    rho_g = 0.0433 * spgr * pres / (zfact * (temp+459.67));
    x = 3.448 + (986.4/(temp+459.67)) + 0.01009*spgr*28.9586;
    y = 2.4 - 0.2*x;
    k = (9.379 + 0.0160 * 28.9586 * spgr) * ((temp+459.67)**1.5);
    k /= 209.2 + 19.26*28.9586*spgr + (temp+459.67);
    mu_g = k * 0.0001 * Math.exp(x * (rho_g**y));
    rho_g *= 62.428;

    gasprops = {gaszfactor:zfact, gasformationvolumefactor:b_g, gasdensity:rho_g, gasviscosity:mu_g};

    return gasprops;
}
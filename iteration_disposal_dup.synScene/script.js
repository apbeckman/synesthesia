var synval = 0;
var bassval = 0;
var midval = 0;
var midhighval = 0;
var hitval = 0;
function update(dt){
    synval = synval + syn_Presence*0.01*spin;
    bassval = bassval + syn_Presence*0.01+ syn_BassPresence*0.03*infl+syn_BassPresence*0.03*rate;
    midval = midval + syn_Presence*0.01+ syn_MidPresence*0.03*infl+syn_MidPresence*0.03*rate;
    midhighval = midhighval + syn_Presence*0.01 + syn_MidHighPresence*0.03;
    hitval = hitval + syn_BassLevel*tunnel_rate*tunnel;
    if (reset_tun > 0.5) hitval = 0;

    setUniform("synvals", synval);
    setUniform("bassvals", bassval);
    setUniform("midvals", midval);
    setUniform("midhighvals", midhighval);
    setUniform("hitvals", hitval);

}
int iFrame = int(FRAMECOUNT);
float time = TIME;
vec2 resolution = RENDERSIZE;

vec3 palette( in float t, in vec3 a, in vec3 b, in vec3 c, in vec3 d )
{
    return a + b*cos( 6.28318*(c*t+d) );
}

  // finalCol *= palette(fract(redAmt + grnAmt + bluAmt), vec3(0.420, 0.500, 0.500), vec3(0.500, 0.500, 0.500), vec3(0.600, 0.250, 1.200), vec3(0.500, 0.450, 0.500));

vec3 rotate(vec3, float);
float distance_field(vec2);
float p_circle(vec2, vec3);
vec2 polar_cart(vec2);
vec2 cart_polar(vec2);
vec3 cwarp(vec3);
vec3 rotate(vec3 c, float a) {
    float ca = cos(a);
    float sa = sin(a);
    return vec3((((c).x) * (ca)) - (((c).y) * (sa)), (((c).y) * (ca)) + (((c).x) * (sa)), (c).z);
}
float distance_field(vec2 p) {
    float dist = max(p_circle(p, vec3(0.0, 0.0, 0.6+fract(p.x*10.0)*distort)), - (max(max(max(max(p_circle(p, vec3(0.0, 0.0, 0.4)), - (max(max(p_circle(p, rotate(vec3(0.13625, -0.04, 0.52221283975), 4.18878666667)), p_circle(p, rotate(vec3(-0.25375, 0.135, 0.464630229322), 4.18878666667))), - (p_circle(p, rotate(vec3(-0.2925, -0.10125, 0.46492270863), 4.18878666667)))))), - (max(max(p_circle(p, rotate(vec3(0.13625, -0.04, 0.52221283975), 2.09439333333)), p_circle(p, rotate(vec3(-0.25375, 0.135, 0.464630229322), 2.09439333333))), - (p_circle(p, rotate(vec3(-0.2925, -0.10125, 0.46492270863), 2.09439333333)))))), - (max(max(p_circle(p, rotate(vec3(0.13625, -0.04, 0.52221283975), 0.0)), p_circle(p, rotate(vec3(-0.25375, 0.135, 0.464630229322), 0.0))), - (p_circle(p, rotate(vec3(-0.2925, -0.10125, 0.46492270863), 0.0)))))), - (max(max(max(p_circle(p, vec3(0.0, 0.18375, 0.295047665641)), p_circle(p, rotate(vec3(0.0, 0.18375, 0.295047665641), 2.09439333333))), p_circle(p, rotate(vec3(0.0, 0.18375, 0.295047665641), 4.18878666667))), - (p_circle(p, vec3(0.0, 0.0, 0.0434))))))));
    float sign = (dist) / (dist);
    return (sin(((dist) * (50.0 + syn_Presence*50.0)) + ((syn_HighTime*0.25) * 5.0))) * (sign);
}
float p_circle(vec2 p, vec3 c) {
    c = cwarp(c);
    return (length(((c).xy) - (p))) - ((c).z);
}
vec2 polar_cart(vec2 p) {
    return (vec2(cos((p).x), sin((p).x))) * ((p).y);
}
vec2 cart_polar(vec2 p) {
    return vec2(atan((p).y, (p).x), length(p));
}
vec4 pattern(vec2 fragCoord) {
    vec2 p = ((((((fragCoord).xy) / ((RENDERSIZE).xy)) * 2.0) - 1.0) * (vec2(1.0, ((RENDERSIZE).y) / ((RENDERSIZE).x)))) * (vec2(1.2, -1.2));
    vec2 h = vec2(0.001, 0.0);
    float d = clamp(((abs(distance_field(p))) / (length((vec2((distance_field((p) + (h))) - (distance_field((p) - (h))), (distance_field((p) + ((h).yx))) - (distance_field((p) - ((h).yx))))) / (2.0 * ((h).x))))) * 400.0, 0.0, 1.0);
    return vec4(d, d, d, 1.0);
}
vec3 cwarp(vec3 c) {
    float t = circle_or_tri;
    return vec3(polar_cart((cart_polar(((c).xy) * (1.0 - (t)))) + (vec2(((t) * 3.14159) * -1.0, 0.0))), ((c).z) * (((t) * 2.5) + 1.2));
}

vec4 renderMain()
{
    if (PASSINDEX == 0.0){
    vec2 uv = _xy.xy / RENDERSIZE;
    vec4 fragColor = vec4(0.0);
    fragColor = 1.0-pattern(_xy);
    if (syn_MediaType>0.5){
        fragColor *= mix(vec4(1.0),_loadUserImage(),media_on);
    }

     

        
        // use texture channel0 for color? why not.
        // vec3 cexp = texture(colornoise, uv * 0.001).xyz * 3.0 + texture(colornoise, uv * 0.01).xyz;//vec3(1.0, 2.0, 4.0);
        // cexp *= 1.4;
        // vec3 col = vec3(pow(v, cexp.x), pow(v, cexp.y), pow(v, cexp.z)) * 2.0;
        // fragColor = vec4(col, 1.0);

        return fragColor;
    }
    else if (PASSINDEX == 1.0){
        vec4 finalCol = vec4(0.0);
        const int iters = 5;
        vec4 last = vec4(1.0);
        vec4 burn = vec4(0.0);
        for(int j = 0; j < iters; j++){
            
            float index = mod((TIME*0.1+syn_BassTime*0.25) + j,iters);
            // float zoom = 0.33+index*0.33;
            float zoom = pow(1.0-index/iters,0.05+power_zoom);

            float divisor = zoom;
            vec2 uv2 = _uv;
            uv2 -= 0.5;
            // uv2 *= 2.0;
            uv2 *= zoom;
            uv2 *= vec2(RENDERSIZE.x/RENDERSIZE.y, 1.0);
            uv2 = _rotate(uv2, j*0.2*syn_BassPresence*twist);
            uv2 /= vec2(RENDERSIZE.x/RENDERSIZE.y, 1.0);

            // uv2 += subber*0.1;
            uv2 += 0.5;
            // uv2 *= 0.5;


            float mixer = (1.0-smoothstep(iters-1, iters, index))*smoothstep(0, 1, index);
            vec4 new = texture(backbuffer, uv2);
            last = new;
            finalCol += new*mixer;
        }

        finalCol /= (iters*0.5);
        if (syn_MediaType > 0.5){
        	finalCol *= 2.0;
        	finalCol = clamp(finalCol, 0.0, 1.0);
        }

        // finalCol = pow(burn,vec4(3.0));
        finalCol += pow(finalCol,vec4(5.0))*syn_HighHits*flashing*100.0;

        finalCol = mix(finalCol, texture(backbuffer, _uv), 1.0-layers_on);

        return finalCol;
    }
}
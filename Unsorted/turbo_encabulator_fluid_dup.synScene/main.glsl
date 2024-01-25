vec3 _grad3(vec3 col1, vec3 col2, vec3 col3, float mixVal){
    mixVal *= 2.0;
    float mix1 = clamp(mixVal,0.0,1.0);
    float mix2 = clamp(mixVal-1.0, 0.0, 1.0);
    return mix(mix(col1, col2, mix1), mix(col2, col3, mix2), step(1.0, mixVal));
}

vec3 _grad4(vec3 col1, vec3 col2, vec3 col3, vec3 col4, float mixVal){
    mixVal *= 3.0;
    float mix1 = clamp(mixVal,0.0,1.0);
    float mix2 = clamp(mixVal-1.0, 0.0, 1.0);
    float mix3 = clamp(mixVal-2.0, 0.0, 1.0);
    vec3 firstTwo = mix(mix(col1, col2, mix1), mix(col2, col3, mix2), step(1.0, mixVal));
    return mix(firstTwo, mix(col3, col4, mix3), step(2.0, mixVal));
}

			//******** BuffA Code Begins ********
vec3 components;
vec3 tx;
vec3 txPlain;

vec2 hash( vec2 p ) // replace this by something better
{
	p = vec2( dot(p,vec2(127.1,311.7)),
			  dot(p,vec2(269.5,183.3)) );

	return -1.0 + 2.0*fract(sin(p)*43758.5453123);
}

// Begin IQ's simplex noise:

// The MIT License
// Copyright © 2013 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
float noise( in vec2 p )
{
    const float K1 = 0.366025404; // (sqrt(3)-1)/2;
    const float K2 = 0.211324865; // (3-sqrt(3))/6;

	vec2 i = floor( p + (p.x+p.y)*K1 );

    vec2 a = p - i + (i.x+i.y)*K2;
    vec2 o = step(a.yx,a.xy);
    vec2 b = a - o + K2;
	vec2 c = a - 1.0 + 2.0*K2;

    vec3 h = max( 0.5-vec3(dot(a,a), dot(b,b), dot(c,c) ), 0.0 );

	vec3 n = h*h*h*h*vec3( dot(a,hash(i+0.0)), dot(b,hash(i+o)), dot(c,hash(i+1.0)));

    return dot( n, vec3(70.0) );

}

// End IQ's simplex noise

//arc flashing functions
float uvrand(vec2 uv)
{
    return fract(sin(dot(uv, vec2(12.9898, 78.233))) * 43758.5453);
}

float arc(vec2 coord)
{
    const float pi = 3.1415926;
    float t = syn_BeatTime * 7.3962;

    vec2 sc = coord;
    float phi = atan(sc.y, sc.x + 1e-6);
    vec2 pc = vec2(fract(phi / (pi * 2) + syn_BassTime*0.01), length(sc));

    vec2 org = vec2(0.5+0.5*sin(syn_BeatTime), 0.5+0.5*sin(syn_BeatTime));
    vec2 wid = vec2(0.5, 0.5);

    for (int i = 0; i < 4; i++)
    {
        if (uvrand(org + t) < 0.04 * i) break;
        wid *= 0.5;
        org += wid * (step(org, pc) * 2 - 1);
    }

    return uvrand(org);
}

float arcGridPattern()
{
    vec4 delta = vec4(-1, -1, 1, 1) * 0.5;
    delta *= (1/RENDERSIZE.y);

    // neightbor four samples
    float c1 = arc(abs(_uvc)*0.4 + delta.xy);
    float c2 = arc(abs(_uvc)*0.4 + delta.zy);
    float c3 = arc(abs(_uvc)*0.4 + delta.xw);
    float c4 = arc(abs(_uvc)*0.4 + delta.zw);

    // roberts cross operator
    float gx = c1 - c4;
    float gy = c2 - c3;
    float g = sqrt(gx * gx + gy * gy);

    return g*4;
}


vec2 pattern(vec2 p)
{
  p = fract(p);
  float r = 10.123;
  float v = 0.0, g = 0.0;
  r = fract(r * 9184.928);
  float cp, d;

  d = p.x;
  g += pow(clamp(1.0 - abs(d), 0.0, 1.0), 1000.0);
  d = p.y;
  g += pow(clamp(1.0 - abs(d), 0.0, 1.0), 1000.0);
  d = p.x - 1.0;
  g += pow(clamp(3.0 - abs(d), 0.0, 1.0), 1000.0);
  d = p.y - 1.0;
  g += pow(clamp(1.0 - abs(d), 0.0, 1.0), 10000.0);

  const int iter = 12;
  for(int i = 0; i < iter; i ++)
  {
    cp = 0.5 + (r - 0.5) * 0.9;
    d = p.x - cp;
    g += pow(clamp(1.0 - abs(d), 0.0, 1.0), 200.0);
    if(d > 0.0) {
      r = fract(r * 4829.013);
      p.x = (p.x - cp) / (1.0 - cp);
      v += 1.0;
    }
    else {
      r = fract(r * 1239.528);
      p.x = p.x / cp;
    }
    p = p.yx;
  }
  v /= float(iter);
  return vec2(g, v);
}

vec2 normz(vec2 x) {
	return x == vec2(0.0, 0.0) ? vec2(0.0, 0.0) : normalize(x);
}

// reverse advection
vec3 advect(sampler2D buff, vec2 ab, vec2 vUv, vec2 step, float sc) {

    vec2 aUv = vUv - ab * sc * step;

    const float _G0 = 0.25; // center weight
    const float _G1 = 0.125; // edge-neighbors
    const float _G2 = 0.0625; // vertex-neighbors

    // 3x3 neighborhood coordinates
    float step_x = step.x;
    float step_y = step.y;
    vec2 n  = vec2(0.0, step_y);
    vec2 ne = vec2(step_x, step_y);
    vec2 e  = vec2(step_x, 0.0);
    vec2 se = vec2(step_x, -step_y);
    vec2 s  = vec2(0.0, -step_y);
    vec2 sw = vec2(-step_x, -step_y);
    vec2 w  = vec2(-step_x, 0.0);
    vec2 nw = vec2(-step_x, step_y);

    vec3 uv =    texture(buff, fract(aUv)).xyz;
    aUv -= _pulse(uv.z, -0.8, 0.2)*ab*0.2*syn_BPMTri2/syn_BPMConfidence;
    vec3 uv_n =  texture(buff, fract(aUv+n)).xyz;
    vec3 uv_e =  texture(buff, fract(aUv+e)).xyz;
    vec3 uv_s =  texture(buff, fract(aUv+s)).xyz;
    vec3 uv_w =  texture(buff, fract(aUv+w)).xyz;
    vec3 uv_nw = texture(buff, fract(aUv+nw)).xyz;
    vec3 uv_sw = texture(buff, fract(aUv+sw)).xyz;
    vec3 uv_ne = texture(buff, fract(aUv+ne)).xyz;
    vec3 uv_se = texture(buff, fract(aUv+se)).xyz;

    return _G0*uv + _G1*(uv_n + uv_e + uv_w + uv_s) + _G2*(uv_nw + uv_sw + uv_ne + uv_se);
}

float pulseMedia = media_pulse_auto*((syn_BassLevel*0.5+syn_Intensity*0.5));
float sizeMod = mix(1.0, 0.5+(sin(_uvc.x*5.0+TIME*0.1)*_uv.y+cos(_uvc.x*4.7-TIME*0.47)+cos(_uvc.x*11.0+TIME*0.89)*(1.0-_uv.y))/6.0, 1.0);
float dist(vec2 p0, vec2 pf){return sqrt((pf.x-p0.x)*(pf.x-p0.x)+(pf.y-p0.y)*(pf.y-p0.y));}

float d = dist(RENDERSIZE.xy*0.5,_xy.xy)*(_mouse.x/RENDERSIZE.x+0.1)*0.005;

float d2 = dist(RENDERSIZE.xy*0.5,_xy.xy)*0.00125;

vec4 renderPassA(sampler2D buff) {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

  const float _K0 = -20.0/6.0; // center weight
  const float _K1 = 4.0/6.0;   // edge-neighbors
  const float _K2 = 1.0/6.0;   // vertex-neighbors
	float cs = curl_scale; // curl scale

  float ls = 0.3-syn_BassPresence*0.295;  // laplacian scale
  float ps = meltyDetail - (cos(syn_BeatTime*0.1) * 0.5);  // laplacian of divergence scale
  float ds = -abs(0.13 - (syn_ToggleOnBeat/50.)); // divergence scale
  const float dp = -0.04; // divergence update scale
  const float pl = 0.3;   // divergence smoothing
  float ad = 4. + syn_Level*2.0 - syn_HighHits*2.0;   // advection distance scale
  ad = mix(ad, 6.0*manual_growth, 1.0-auto_growth);

  const float pwr = 1.0;  // power when deriving rotation angle from curl
  const float amp = 1.0;  // self-amplification
  float upd = 0.7 - (syn_BassHits/3.)*0.6;  // update smoothing
  upd = mix(upd, 0.9-manual_growth*0.6, 1.0-auto_growth);
  float sq2 = 0.6 + scatter;  // diagonal weight



	vec2 vUv = _uv;
  vec2 texel = sizeMod*size_modifier / RENDERSIZE.xy;
	vUv -= _uvc*0.0005*zooming;
  vUv += _uvc*d2*0.0005*fisheye;
  // 3x3 neighborhood coordinates
  float step_x = texel.x;
  float step_y = texel.y;
  vec2 n  = vec2(0.0, step_y);
  vec2 ne = vec2(step_x, step_y);
  vec2 e  = vec2(step_x, 0.0);
  vec2 se = vec2(step_x, -step_y);
  vec2 s  = vec2(0.0, -step_y);
  vec2 sw = vec2(-step_x, -step_y);
  vec2 w  = vec2(-step_x, 0.0);
  vec2 nw = vec2(-step_x, step_y);

	// uv.x and uv.y are the x and y components, uv.z is divergence
  vec3 uv =    texture(buff, fract(vUv)).xyz;
  vec3 uv_n =  texture(buff, fract(vUv+n)).xyz;
  vec3 uv_e =  texture(buff, fract(vUv+e)).xyz;
  vec3 uv_s =  texture(buff, fract(vUv+s)).xyz;
  vec3 uv_w =  texture(buff, fract(vUv+w)).xyz;
  vec3 uv_nw = texture(buff, fract(vUv+nw)).xyz;
  vec3 uv_sw = texture(buff, fract(vUv+sw)).xyz;
  vec3 uv_ne = texture(buff, fract(vUv+ne)).xyz;
  vec3 uv_se = texture(buff, fract(vUv+se)).xyz;
  // float pat = pattern(_toPolar((_uv-0.5)*vec2(RENDERSIZE.x/RENDERSIZE.y,1.0))*vec2(1.0, 1.0)+vec2(-syn_BassTime*0.005, 0.0)).r;
  // pat = clamp(step(1.4, pat), 0.0, 1.0)*_pulse(length(_uvc), 0.5, 0.35);
  vec3 blurUv = (uv_n+uv_s+uv_w+uv_e)/4.0;

  float pat = arcGridPattern();

	vec3 arcFlash = vec3(1.0)*pat*syn_OnBeat*disrupt_field;

    // laplacian of all components
    vec3 lapl  = _K0*uv + _K1*(uv_n + uv_e + uv_w + uv_s) + _K2*(uv_nw + uv_sw + uv_ne + uv_se);
    float sp = ps * lapl.z;

    // calculate curl
    // vectors point clockwise about the center point
    float curl = uv_n.x - uv_s.x - uv_e.y + uv_w.y + sq2 * (uv_nw.x + uv_nw.y + uv_ne.x - uv_ne.y + uv_sw.y - uv_sw.x - uv_se.y - uv_se.x);

    // compute angle of rotation from curl
    float sc = cs * sign(curl) * pow(abs(curl), pwr);


    // calculate divergence
    // vectors point inwards towards the center point
    float div  = uv_s.y - uv_n.y - uv_e.x + uv_w.x + sq2 * (uv_nw.x - uv_nw.y - uv_ne.x - uv_ne.y + uv_sw.x + uv_sw.y + uv_se.y - uv_se.x);
    float sd = uv.z + dp * div + pl * lapl.z;

    vec2 norm = normz(uv.xy);
    vec3 ab = advect(buff, vec2(uv.x, uv.y), vUv, texel, ad);

    // temp values for the update rule
    float ta = amp * ab.x + ls * lapl.x + norm.x * sp + uv.x * ds * sd;
    float tb = amp * ab.y + ls * lapl.y + norm.y * sp + uv.y * ds * sd;

    // rotate
    float a = ta * cos(sc) - tb * sin(sc);
    float b = ta * sin(sc) + tb * cos(sc);

    uv = mix(uv, blurUv, syn_BPMTri2/syn_BPMConfidence);
    float addSpot = 1/(0.5+distance(_uvc, vec2(0.5*sin(TIME), 0.5*cos(TIME))));
    // uv.z += (-0.5+(syn_BPMTri2/syn_BPMConfidence))*0.05*addSpot*_pulse(uv.z, 0.0, 0.1);
    vec3 abd = upd*uv + (1.0-upd)*vec3(a,b,sd);

		//add arc flashing
		abd.xy += (normz(arcFlash.xy))*mix(-1.0, 1.0, syn_OnBeat)*normalize(_uvc);

		float mediaLum = 0.5-(dot(txPlain.rgb, vec3(1.0))/3.0);

		abd.xy -= 1.0*mediaLum*pulseMedia;

    abd.z -= mediaLum*media_hint*media_hint*media_hint*media_hint*0.5;

		//Tracers
		if(Tracers_on > 0.5)
		{
			vec2 d = fragCoord.xy;
			float maxTracers = 8.;
			vec2 flowVec = vec2(0.);
			for(int i = 0; i < maxTracers; i++)
			{
				flowVec = vec2(abs(sin(syn_Time*2*PI*0.004*(syn_Presence*sin(i))*tracer_intensity)),abs(cos(syn_HighTime*2*PI*0.008*(syn_HighPresence*cos(i))+i*tracer_intensity)));
				flowVec += tracer_control;
				d  = fragCoord.xy - smoothstep(-0.1, 1.1, flowVec)*RENDERSIZE.xy;
				float m = exp(-length(d) / (tracer_size + i - syn_BassHits));
				abd.xy += m * normz(d);
			}
		}

    // initialize with noise
    if(FRAMECOUNT<3 || reset > 0.5) {
		vec3 rnd = vec3(noise(7.0 * vUv + 1.1 * syn_HighLevel), noise(9.0 * vUv + 2.2 * syn_BassLevel), noise(8.0 * vUv + 3.3 * syn_MidLevel));
			fragColor = vec4(rnd, 1);
    } else {
      abd.z = clamp(abd.z, -1.0, 1.0);
      abd.xy = clamp(length(abd.xy) > 1.0 ? normz(abd.xy) : abd.xy, -1.0, 1.0);
      abd.xy = abd.xy;
      fragColor = vec4(abd, 1.0);
    }

	return fragColor;
 }

// Visualization of the system in Buffer A
vec4 renderMainImage(sampler2D buff) {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
    float freq = 5.0;
    // vec2 warp = 0.5000*cos( uv.xy*1.0*freq + vec2(0.0,1.0) + TIME*0.2 ) +
    //             0.2500*cos( uv.yx*2.3*freq + vec2(1.0,2.0) + TIME*0.2 ) +
    //             0.1250*cos( uv.xy*4.1*freq + vec2(5.0,3.0) + TIME*0.2 ) +
    //             0.0625*cos( uv.yx*7.9*freq + vec2(3.0,4.0) + TIME*0.2 );

    // uv += warp*0.01;
    vec3 c = texture(buff, uv).xyz;
    vec3 norm = normalize(c);
		//Color assignment
		//TODO: fine tune colors
		float color1 = _pulse(color_palette, 0.0, 1.0);
		float color2 = _pulse(color_palette, 1.0, 1.0);
		float color3 = _pulse(color_palette, 2.0, 1.0);

    vec3 rbcol =  0.2 - 0.75*cross(norm.xyz, vec3(color1, -color2, color3));
		// vec3 rbcol = 0.4 + color_tuning * cross(norm.xyz, vec3(0.5, -0.4, 0.5));
    // rbcol = vec3(0.0);
    // vec3 texCol = _loadMedia().rgb;
    norm.xy = _rotate(norm.xy, syn_HighTime*0.4*color_rotate);
    // float lum = pow(dot(texCol, vec3(1.0))/3.0,2.0)*2.0;
    // lum = clamp(lum, 0.0, 1.0);
    rbcol += norm.x*vec3(0.1,0.6,0.7)*(1.0+0.3*syn_HighHits)*0.5;
    rbcol += norm.y*vec3(0.9,0.6,0.0)*0.5;
    rbcol += vec3(pow(dot(norm,vec3(light_pos, 0.0)),3.0))*0.75;
    // rbcol *= mix(1.0, pow(clamp(c.b, 0.0, 1.0), 0.15), dark_valleys);
    // rbcol += texCol;
    // col += texCol*texCol*0.25;

    // col += f*pow(threeColMix(vec3(0.0), vec3(1.0,0.0,0.0), vec3(0.0,0.0,1.0), pow(f,10.0)),vec3(2.0));
    // col += _pulse(f, 1.0, 0.01)*syn_HighHits*noise(vec3(_uvc*15.0,TIME*2.0))*1.5;
    // col += 2.0*f*pow(_palette(diff, vec3(0.500, 0.500, 0.500), vec3(0.500, 0.500, 0.500), vec3(0.500, 0.250, 0.750), vec3(0.500, 0.500, 0.500)),vec3(2.0));


		//attempts at color tuning
    // rbcol = _rgb2hsv(rbcol);
		// rbcol.x -= color3/2;
		// fragColor = vec4(_hsv2rgb(rbcol), 1.);

		fragColor = vec4(rbcol, 1.);
		fragColor = _gamma(fragColor, 1.1);
		fragColor = clamp(fragColor, 0.0, 1.0);
		//mix in media for final image

    vec3 texCol = vec3(0.0);

    if (media_multiply > 0.5){
        texCol = fragColor.rgb*tx*1.2;
    } else {
        texCol = 2.0*tx*dot(fragColor.rgb, vec3(1.0))/3.0;
    }

    fragColor.rgb = mix(fragColor.rgb, texCol, media_color_mix);
	return fragColor;
 }


vec4 renderMain(){
  components = texture(BuffA, _uv).xyz;
  tx = _loadMedia(vec2(components.z*media_refraction*media_refraction)).xyz;
  txPlain = _loadMedia().xyz;

	if(PASSINDEX == 0){
		return renderPassA(BuffB);
	}
	if(PASSINDEX == 1){
		if (run_twice_fps > 0.5){
			return renderPassA(BuffA);
		} else {
			return texture(BuffA, _uv);
		}
	}
	if(PASSINDEX == 2){
		return renderMainImage(BuffA);
	}
}

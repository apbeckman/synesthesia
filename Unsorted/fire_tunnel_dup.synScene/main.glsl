vec3 filter_() {
  vec2 delta = 1. / RENDERSIZE;
  
  vec3 val = texture(tunnelPass, _uv).xyz;


  vec3 l = texture(tunnelPass, _uv + vec2(0., delta.y)).xyz;
  vec3 r = texture(tunnelPass, _uv - vec2(0., delta.y)).xyz;
  vec3 u = texture(tunnelPass, _uv + vec2(delta.x, 0.)).xyz;
  vec3 d = texture(tunnelPass, _uv - vec2(delta.x, 0.)).xyz;

  vec3 n = vec3(_rand(_uvc+fract(smoothTimeC))) - 0.5;
  
  vec3 bloom = max(val, max(max(l, r), max( u, d)));
  bloom = bloom  + l + r + u + d;
  bloom /= 5.; // orlando;
  return bloom + n/9.;

}

// Created by Stephane Cuillerdier - Aiekick/2015
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// Tuned via XShade (http://www.funparadigm.com/xshade/)

/* 
	variation more cloudy off Another Cloudy Tunnel 3 : 
		https://www.shadertoy.com/view/ltBSzV

	the cloudy famous tech come from the shader of duke : https://www.shadertoy.com/view/MljXDw
        Himself a Port of a demo by Las => http://www.pouet.net/topic.php?which=7920&page=29&x=14&y=9
*/

float t;
float cosPath(vec3 p, vec3 dec){return dec.x * cos(p.z * dec.y + dec.z);}
float sinPath(vec3 p, vec3 dec){return dec.x * sin(p.z * dec.y + dec.z);}
#define Tunnel(p,pos,c,s) p.xy-pos-vec2(cosPath(p, c)*0.5,sinPath(p, s)*0.5)

vec3 getHotColor(float Temp)
{
	vec3 col = vec3(255.);
	col.x = 56100000. * pow(Temp,(-3. / 2.)) + 148.;
   	col.y = 100.04 * log(Temp) - 623.6;
   	if (Temp > 6500.) col.y = 35200000. * pow(Temp,(-3. / 2.)) + 184.;
   	col.z = 194.18 * log(Temp) - 1448.6;
   	col = clamp(col, 0., 255.)/255.;
	if (Temp < 1000.) col *= Temp/1000.;
   	return col;
}

float pn( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);
	f = f*f*(3.0-2.0*f);
	vec2 uv = (p.xy+vec2(37.0,17.0)*p.z) + f.xy;
	vec2 rg = textureLod(image30, (uv+ 0.5)/256.0, 0.0 ).yx;
	return -1.0+2.4*mix( rg.y, rg.x, f.z );
}
float planesMod = pow(planes_only, 10.0);
float df(vec3 p)
{
	float pnNoise = pn(p*.26)*1.98 + pn(p*.26)*.62 + pn(p*1.17)*.39;
	float path = sinPath(p ,vec3(5.704,0.3828,0.16));
    float bt = min(p.y, -p.y) + 12.;
	float df;
    vec4 vec, var;
	for (float i=0.;i<5.;i++)
	{
		var = vec4(0.35,2.5,5.8,0) + vec4(0.5,0.06,0,5) * i;
		vec.xy = Tunnel(p, vec2(path, 0.), var.xyz, var.zxy);
		df = var.w - min(length(vec.xy), length(vec.zw));
		vec.zw = mix(vec.zw, vec.xy, planesMod);	
	}
	return min(df, bt) + pnNoise*noise;
}

vec3 cam(vec2 uv, vec3 ro, vec3 cu, vec3 cv, float fov)
{
	vec3 rov = normalize(cv-ro);
    vec3 u = normalize(cross(cu, rov));
    vec3 v = normalize(cross(rov, u));
    vec3 rd = normalize(rov + fov*u*uv.x + fov*v*uv.y);
    return rd;
}

vec4 march(vec4 f, vec3 ro, vec3 rd, float st)
{
	float s = 1., h = .25*constriction, td = 0.0, d=1.,dl=0.2, w;
	vec3 p = ro;
	for(float i=0.;i<240.0;i++)
	{      
		p.xy = mix(p.xy, _rotate(p.xy, p.z), twist);
		if(s<0.01||d>60.||td>.95) break;
        s = df(p) * .09 * i/60.;
		if (s < h)
		{
			w = (1.-td) * (h-s) * i/60.;
			f.rgb += texture(syn_UserImage, p.xz*0.01+p.y*0.01).rgb*0.02*(syn_HighPresence+syn_HighHits*2.0)*_pulse(fract(p.y*0.1-(smoothTimeB*0.5)), 0.5, 0.5);
			f.rgb += getHotColor(td*i*70.)*mix(w, (0.05-w), syn_BassPresence*(syn_BassLevel+syn_BassHits)*_pulse(fract(p.z*0.1-(smoothTimeB*0.5)), 0.5, 0.5));
			td += w;
		}
		dl += 1.01 - exp(-0.001 * log(d));	
		td += 0.01;
		s = max(s, st);
		d+=s;
		p = ro+rd*d;	
   	}
	f.rgb = mix( f.rgb, vec3(0), 1.0 - exp( -0.002*d*d) )/dl; // fog
	return f;
}

vec4 renderMainImage() {
	vec4 f = vec4(0.0);
	vec2 g = _xy;

		t = smoothTime+const_speed;
		f = vec4(syn_MidPresence*0.25,0.15,0.32,1);
	  vec2 q = g/RENDERSIZE.xy;
	  vec3 ro = vec3(cos(t*.3), sin(t*.2),t )*vec3(8.5, 8.5, 5.)*0.4;
		vec3 rd = cam((2.*g-RENDERSIZE.xy)/RENDERSIZE.y, ro, vec3(0,1,0), ro + vec3(0,0,1), 3.5*fov);

    rd.xy = _rotate(rd.xy, -spin*PI);

		f = march(f, ro, rd, sharpness+0.01);
	  f.rgb *= 0.5 + 0.5*pow( 16.0*q.x*q.y*(1.0-q.x)*(1.0-q.y), 0.25 ); // iq vignette

	return f; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
	if (PASSINDEX == 1){
		return vec4(filter_(), 1.0);
	}
}

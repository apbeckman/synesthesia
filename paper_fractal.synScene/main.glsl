


// Code by Flopine
// Thanks to wsmind, leon, XT95, lsdlive, lamogui, Coyhot, Alkama and YX for teaching me
// Thanks LJ for giving me the love of shadercoding :3

// Thanks to the Cookie Collective, which build a cozy and safe environment for me 
// and other to sprout :)  https://twitter.com/CookieDemoparty

float hash21 (vec2 x)
{return fract(sin(dot(x,vec2(54.4,62.1)))*457.5);}

mat2 rot (float a)
{return mat2 (cos(a),sin(a*Slider1),-sin(a*Slider2),cos(a));}

void mo (inout vec2 p, vec2 d)
{
    p = abs(p)-d;
    if (p.y>p.x) p = p.yx;
}

float plane (vec3 p, vec3 n)
{return dot(p, normalize(n));}

float cut_ps (vec3 p, float s)
{
    p *= s;
    mo(p.xy,vec2(1.));
    mo(p.yz,vec2(0.6*Slider8+Slider8B));
    mo(p.xz, vec2(0.1*Slider8+Slider8B));
    return plane(p,vec3(1.,1.,4.))/(s);
}

float prim1 (vec3 p, float s)
{
    p.xz *= rot(TIME*Slider3+Slider5);
    return cut_ps(p,s);
}

float fractal(vec3 p)
{
    float size = 1.;
    float d = prim1(p,size-Slider7);
    for (int i=1; i<5; i++)
    {
        float ratio = float(i)/5.;         
        p.yz *= rot(TIME*ratio*Slider4+Slider6);
		size -= 0.2;
        d = min(d, prim1(p,size+Slider7));
    }
    return d;
}

float g1 = 0.;
float SDF (vec3 p)
{
    float noise = IMG_NORM_PIXEL(iChannel0,mod((p.xy*0.1)+TIME,1.0)).r;
    float sphe = length(p)-(.8+noise);
    g1 += 0.1/(0.1+sphe*sphe);
    return max(-length(p+vec3(0.,0.,4.5))+.8,min(sphe,fractal(p)));
}

vec4 renderMain() { 
 	vec4 out_FragColor = vec4(0.0);




    vec2 uv = (2.*_xy.xy-RENDERSIZE.xy)/RENDERSIZE.y;
    float dither = hash21(uv);
    
    vec3 ro = vec3(0.001+10.-Sliderx,0.001+5.-Slidery,-4.5-Sliderz),
        rd = normalize(vec3(uv,0.8)),
        p = ro,
        col = vec3(0.1);
    
    float shad = 0.;
    bool hit = false;
    for (float i=0.; i<64.; i++)
    {
        float d = SDF(p);
        if (d<0.001)
        {
            hit = true;
            shad = i/64.;
            break;
        }
        d *= 0.8+dither*0.1;
        p += d*rd;
    }
    if (hit)
    {
        col = vec3(1.-shad);
		col += g1*vec3(0.15,0.,0.1);
    }
    
    out_FragColor = vec4(col,1.0);

return out_FragColor; 
 } 


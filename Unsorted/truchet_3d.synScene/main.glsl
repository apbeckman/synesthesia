vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


float N21(vec2 p){
    p+=12.245;
    p=fract(p*vec2(1.432342,515.5621));
    p+=dot(p,p+341.3535);
    return fract(p.x*p.y*.12);
    
}


float sdTorus( vec3 p, vec2 t )
{
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  float r=length(q)-t.y;
  return min(r,length(vec3(abs(p.x),p.y,abs(p.z))-vec3(.17,0,.18))-.06);
}


float sdBox( vec3 p, vec3 b )
{
  vec3 q = abs(p) - b;
  
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float crossTorus(vec3 p,vec3 id){
    float w=.035;
    
    float t1=sdTorus(vec3(abs(p.x)-.5,p.z,abs(p.y)-.5),vec2(.25,w));
    //t1=min(t1,length(p-vec3(.25,-.25,-.25))-.2);
       


    float t2=sdTorus(vec3(abs(p.x)-.5,p.z,p.y),vec2(.25,w));
    t2=min(t2,sdTorus(vec3(p.x,p.z,abs(p.y)-.5),vec2(.25,w)));
   


    float r=N21(id.xy+id.z+floor(TIME*0.1));
    float rn=N21(id.xy+id.z+floor(TIME*0.1)+1.);
    float t=r>.5?t1:t2;
    float tn=rn>.5?t1:t2;
    
    return mix(t,tn,smoothstep(.8,.9,fract(smoothTimeC/3.)));

    
}

float cube(vec3 p, vec3 id){
    
    float c=cos(1.57);
    float s=sin(1.57);
    mat2 m=mat2(c,s,-s,c);
    vec2 a=p.yz*m;
    vec2 b=p.xz*m;
    float t=min(crossTorus(p,id),crossTorus(vec3(p.x,a.x,a.y),id));
    t=min(crossTorus(vec3(a.x,p.y,b.y),id),t);
    return max(sdBox(p,vec3(.52)),t);
}

vec3 transform(inout vec3 p){
    float c=cos(sin(TIME/10.));
    float s=sin(sin(TIME/10.));
    
    mat2 m=mat2(c,s,-s,c);
    p.xy*=m;
    vec3 op=p;
    vec3 r=vec3(1.);
    p=fract(p)-.5;
    return floor(op-p); 
}

float dist(vec3 p){
    
    vec3 id=transform(p);   
    float cu=cube(p,id);
       
    return cu;
}

vec3 normal(vec3 p){
    vec2 e=vec2(.001,0.0);
    float d=dist(p);
    vec3 n=vec3(
        d-dist(p-e.xyy),
        d-dist(p-e.yxy),
        d-dist(p-e.yyx));
    return normalize(n);
}



float rayMarch(vec3 ro, vec3 rd){
    float d=0.;
    float td=0.;
    float i=0.;;
    for(;i<50.;i++){
        d=dist(ro+td*rd);
        td+=d;
        if(td>50. || d<.001)
            break;
    }
    if(i>50.)
        return 50.;
    return float(td);    
}

float softshadow( in vec3 ro, in vec3 rd, float k )
{
    float res = 1.0;

    for( float t=.001; t<100.; )
    {
        float h = dist(ro + rd*t);
        if( h<0.001 )
            return 0.0;
        res = min( res, k*h/t );
        t += h;
    }
    return res;
}

vec3 shade(vec3 p,vec3 rd){
    vec3 n=normal(p);
    
    
    vec3 id=transform(p); 
    
    vec3 col=vec3(.15,.55,.80);
    //WIP color mix
    vec3 col1=normalize(.5+.5*sin(smoothTimeB/2.+col*(id.x+id.y*10.)));
    vec3 col2=normalize(.5+.5*sin(smoothTimeB/2.+col*((id.x+1.)+id.y*10.)));
    vec3 col3=normalize(.5+.5*sin(smoothTimeB/2.+col*((id.x-1.)+id.y*10.)));
    vec3 col4=normalize(.5+.5*sin(smoothTimeB/2.+col*((id.x)+(id.y-1.)*10.)));
    vec3 col5=normalize(.5+.5*sin(smoothTimeB/2.+col*((id.x)+(id.y+1.)*10.)));
    col1+=col3*smoothstep(.75,1.,length(p.x-vec2(.25)));
    col1+=col2*smoothstep(.75,1.,length(p.x+vec2(.25)));
    col1+=col4*smoothstep(.75,1.,length(p.y-vec2(.25)));
    col1+=col5*smoothstep(.75,1.,length(p.y+vec2(.25)));
    
    vec3 a=col1/1.5*vec3(abs(dot(n,normalize(vec3(3.,2.,13.)-p))));
    return a;//*softshadow(p+n*.1,normalize(-p+vec3(0,0.,0.-TIME/10.)),2.);
}






void mainImage0( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (fragCoord-.5*RENDERSIZE.xy)/RENDERSIZE.y;

    float t=TIME/10.;
    vec3 camera=vec3(0.,0.,6.-t);
    vec3 lookAt=vec3(sin(TIME/3.),cos(TIME/3.),-t);
    float zoom=.8;
    vec3 f=normalize(camera-lookAt);
    vec3 r=cross(vec3(0,1.,0),f);
    vec3 u=cross(f,r);
    vec3 c=camera-f*zoom;
    vec3 i=c+uv.x*r+uv.y*u;
    vec3 ray=normalize(i-camera);
    
    float d=rayMarch(c,ray);
    
    
    vec3 p=c+ray*d;
    vec3 n=normal(p);
    vec3 col=shade(p,ray);
    /*vec3 rf=reflect(n,ray);
    float df=rayMarch(p+n*.01,rf);
    float fresnel=clamp(pow(1.-dot(ray,n),5.),.2,.8);
    col+=shade(p+rf*df,rf)*.3*fresnel;
      */  
    fragColor = vec4(col*max((1.-d*.15),.05),1.);
}

//from https://shadertoyunofficial.wordpress.com/2021/03/09/advanced-tricks/ 
vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 U = _xy;

    mainImage0(O,U);
    if ( fwidth(length(O)) > .05 ) {  // difference threshold between neighbor pixels
        vec4 o;
        for (int k=0; k < 9; k+= k==3?2:1 )
          { mainImage0(o,U+vec2(k%3-1,k/3-1)/3.); O += o; }
        O /= 9.;
      //O.r++;                        // uncomment to see where the oversampling occurs
    }
	return O; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}
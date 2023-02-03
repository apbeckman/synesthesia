

			//******** BuffA Code Begins ********

// created by florian berger (flockaroo) - 2016
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

// record video histroy in Xnum*Ynum grid

#define Xnum 10
#define Ynum 10

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv0=fragCoord/RENDERSIZE.xy;
    int fr=int(uv0.x*float(Xnum))+int(uv0.y*float(Xnum))*Ynum;
    if(fr!=int(mod(float(FRAMECOUNT),float(Xnum*Ynum)))) 
        fragColor = texture(BuffA,uv0);
    else
        fragColor = texture(syn_UserImage,fract(uv0*vec2(Xnum,Ynum)));
    if(FRAMECOUNT<=5)
        fragColor = texture(syn_UserImage,fract(uv0*vec2(Xnum,Ynum)));
	return fragColor; 
 } 


// created by florian berger (flockaroo) - 2016
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

// flip pixels randomly from old frame to new frame

// tiles of the video history (as in BufA)
#define Xnum 10
#define Ynum 10

// tile size in pixels (proportional sqrt(Res), otherwise pixels get too small in preview)
#define TileSize (10.*sqrt(RENDERSIZE.y/1080.))
// take every DFrame'th frame
#define DFrame 30


#define Res1 RENDERSIZE.xy

float rectMask(float b, float w, vec2 uv)
{
	vec4 e=smoothstep(vec4(-b-1.5*w),vec4(-b+.5*w),vec4(uv,vec2(1)-uv));
    return e.x*e.y*e.z*e.w;
}

vec2 frameUV(int frame, vec2 uv)
{
    frame = int(mod(float(frame+Xnum*Ynum),float(Xnum*Ynum)));
    return (uv+vec2(mod(float(frame),float(Xnum)),frame/Xnum))/vec2(Xnum,Ynum);
}

vec4 getColor(float frame, vec2 uv)
{
    vec4 c1=texture(BuffA,frameUV(int(floor(frame)),uv));
    return c1;
    vec4 c2=texture(BuffA,frameUV(int(ceil(frame)),uv));
    return mix(c1,c2,fract(frame));
}

float getFrameSpec(vec2 duv, vec2 fragCoord)
{
    vec2 uvS=.85*duv*pow(abs(.85*duv),vec2(20.))*530.;
    vec3 n=normalize(vec3(uvS,1.1));
    vec3 v=normalize(vec3((fragCoord-RENDERSIZE.xy*.5)/RENDERSIZE.x*.25,-1));
    vec3 light=vec3(3.5,2.3,-2.6);
    vec3 halfVec=normalize(normalize(-v)+normalize(light));
    return pow(clamp(dot(n,halfVec),0.,1.),10.);
}

float getVign(vec2 fragCoord)
{
	float rs=length(fragCoord-RENDERSIZE.xy*.5)/RENDERSIZE.x;
    return 1.-rs*rs*rs;
}

float circle(vec2 uv, float r)
{
    float l=length(uv-.5);
    return 1.-smoothstep(r-.05,r+.05,l);
}

float circleMask(vec2 uv,float y, float i1, float i2)
{
    float r1=.45;
    float r2=.45;
    r1*=i1;
    r2*=i2;
    return circle(uv+vec2(0,y-1.),r1)+circle(uv+vec2(0,y),r2);
}

    
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    int actFrame=(int(FRAMECOUNT)/DFrame)*DFrame;
    int prevFrame=(int(FRAMECOUNT)/DFrame-1)*DFrame;
    vec4 rand = texture(image30,(floor(fragCoord/TileSize+float(FRAMECOUNT/DFrame)*1.*13.)+.5)/Res1);
    
    vec2 uvQ = (floor(fragCoord/TileSize)*TileSize)/RENDERSIZE.xy;
    vec2 uv = fragCoord/RENDERSIZE.xy;
    vec2 duv = (uv-uvQ)*RENDERSIZE.xy/TileSize;
    
    vec4 c1=getColor(float(actFrame),uvQ);
    vec4 c2=getColor(float(prevFrame),uvQ);
    
    // let the pixels flip at a randomly
    float y=-rand.x*2.+3.*float(FRAMECOUNT-actFrame)/float(DFrame);
    y=clamp(y,0.,1.);
    y*=y;
    float r = fract(rand.y+float(actFrame/DFrame)*.25);
    float thr=1.-duv.y;
    
    // some tests flipping in different directions
    //if (r>.25) thr=1.-duv.x;
    //if (r>.50) thr=duv.x;
    //if (r>.75) thr=duv.y;
    
    // some tests playing with different circle sizes prop to color intensity
    //float i1=dot(vec3(.333),c1.xyz);
    //float i2=dot(vec3(.333),c2.xyz);    
    float i1=1.;
    float i2=1.;
    
	fragColor = mix(c1,c2,smoothstep(y-.1,y+.1,thr));
    fragColor = mix(fragColor,vec4(.2,.3,.4,1),.15-.15*circleMask(duv,y,i1,i2));
    
    // some rectangular ambient around every macro-pixel
    fragColor *= .5+.5*rectMask(.2*(dot(fragColor.xyz,vec3(.333))),.7,duv);
    
    float spec = 0.0
        //+getFrameSpec(duv,fragCoord)
        //+getFrameSpec(fract(duv+vec2(0,y)),fragCoord)
        //+getCricleSpec(fract(duv+vec2(0,y)),fragCoord)
        +clamp((circleMask(duv-.02,y,i1,i2)-circleMask(duv+.02,y,i1,i2)),-.4,1.)
        ;
    
    fragColor.xyz += .5*spec;
    fragColor *= 1.2*getVign(fragCoord);
    //fragColor = texture(BuffA,uv);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderMainImage();
	}
}
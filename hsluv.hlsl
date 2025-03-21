/*
HSLUV-HLSL v1.0
HSLUV is a human-friendly alternative to HSL. ( http://www.hsluv.org ) HLSL port by Roman Ilyin ( https://github.com/romanilyin )
*/

float3 hsluv_intersectLineLine(float3 line1x, float3 line1y, float3 line2x, float3 line2y) {
    return (line1y - line2y) / (line2x - line1x);
}

float3 hsluv_distanceFromPole(float3 pointx,float3 pointy) {
    return sqrt(pointx*pointx + pointy*pointy);
}

float3 hsluv_lengthOfRayUntilIntersect(float theta, float3 x, float3 y) {
    float3 len = y / (sin(theta) - x * cos(theta));
    if (len.r < 0.0) {len.r=1000.0;}
    if (len.g < 0.0) {len.g=1000.0;}
    if (len.b < 0.0) {len.b=1000.0;}
    return len;
}

float hsluv_maxSafeChromaForL(float L){
    float3x3 m2 = float3x3(
         3.2409699419045214  ,-0.96924363628087983 , 0.055630079696993609,
        -1.5373831775700935  , 1.8759675015077207  ,-0.20397695888897657 ,
        -0.49861076029300328 , 0.041555057407175613, 1.0569715142428786  
    );
    float sub0 = L + 16.0;
    float sub1 = sub0 * sub0 * sub0 * .000000641;
    float sub2 = sub1 > 0.0088564516790356308 ? sub1 : L / 903.2962962962963;

    float3 top1   = (284517.0 * m2._m00_m10_m20 - 94839.0  * m2._m02_m12_m22) * sub2;
    float3 bottom = (632260.0 * m2._m02_m12_m22 - 126452.0 * m2._m01_m11_m21) * sub2;
    float3 top2   = (838422.0 * m2._m02_m12_m22 + 769860.0 * m2._m01_m11_m21 + 731718.0 * m2._m00_m10_m20) * L * sub2;

    float3 bounds0x = top1 / bottom;
    float3 bounds0y = top2 / bottom;

    float3 bounds1x =              top1 / (bottom+126452.0);
    float3 bounds1y = (top2-769860.0*L) / (bottom+126452.0);

    float3 xs0 = hsluv_intersectLineLine(bounds0x, bounds0y, -1.0 / bounds0x, float3(0.0, 0.0, 0.0));
    float3 xs1 = hsluv_intersectLineLine(bounds1x, bounds1y, -1.0 / bounds1x, float3(0.0, 0.0, 0.0));

    float3 lengths0 = hsluv_distanceFromPole( xs0, bounds0y + xs0 * bounds0x );
    float3 lengths1 = hsluv_distanceFromPole( xs1, bounds1y + xs1 * bounds1x );

    return  min(lengths0.r,
            min(lengths1.r,
            min(lengths0.g,
            min(lengths1.g,
            min(lengths0.b,
                lengths1.b)))));
}

float hsluv_maxChromaForLH(float L, float H) {

    float hrad = radians(H);

    float3x3 m2 = float3x3(
         3.2409699419045214  ,-0.96924363628087983 , 0.055630079696993609,
        -1.5373831775700935  , 1.8759675015077207  ,-0.20397695888897657 ,
        -0.49861076029300328 , 0.041555057407175613, 1.0569715142428786  
    );
    float sub1 = pow(L + 16.0, 3.0) / 1560896.0;
    float sub2 = sub1 > 0.0088564516790356308 ? sub1 : L / 903.2962962962963;

    float3 top1   = (284517.0 * m2._m00_m10_m20 - 94839.0  * m2._m02_m12_m22) * sub2;
    float3 bottom = (632260.0 * m2._m02_m12_m22 - 126452.0 * m2._m01_m11_m21) * sub2;
    float3 top2   = (838422.0 * m2._m02_m12_m22 + 769860.0 * m2._m01_m11_m21 + 731718.0 * m2._m00_m10_m20) * L * sub2;

    float3 bound0x = top1 / bottom;
    float3 bound0y = top2 / bottom;

    float3 bound1x =              top1 / (bottom+126452.0);
    float3 bound1y = (top2-769860.0*L) / (bottom+126452.0);

    float3 lengths0 = hsluv_lengthOfRayUntilIntersect(hrad, bound0x, bound0y );
    float3 lengths1 = hsluv_lengthOfRayUntilIntersect(hrad, bound1x, bound1y );

    return  min(lengths0.r,
            min(lengths1.r,
            min(lengths0.g,
            min(lengths1.g,
            min(lengths0.b,
                lengths1.b)))));
}

float hsluv_fromLinear(float c) {
    return c <= 0.0031308 ? 12.92 * c : 1.055 * pow(c, 1.0 / 2.4) - 0.055;
}
float3 hsluv_fromLinear(float3 c) {
    return float3( hsluv_fromLinear(c.r), hsluv_fromLinear(c.g), hsluv_fromLinear(c.b) );
}

float hsluv_toLinear(float c) {
    return c > 0.04045 ? pow((c + 0.055) / (1.0 + 0.055), 2.4) : c / 12.92;
}

float3 hsluv_toLinear(float3 c) {
    return float3( hsluv_toLinear(c.r), hsluv_toLinear(c.g), hsluv_toLinear(c.b) );
}

float hsluv_yToL(float Y){
    return Y <= 0.0088564516790356308 ? Y * 903.2962962962963 : 116.0 * pow(Y, 1.0 / 3.0) - 16.0;
}

float hsluv_lToY(float L) {
    return L <= 8.0 ? L / 903.2962962962963 : pow((L + 16.0) / 116.0, 3.0);
}

float3 xyzToRgb(float3 tuple) {
    const float3x3 m = float3x3( 
        3.2409699419045214  ,-1.5373831775700935 ,-0.49861076029300328 ,
       -0.96924363628087983 , 1.8759675015077207 , 0.041555057407175613,
        0.055630079696993609,-0.20397695888897657, 1.0569715142428786  );
    
    return hsluv_fromLinear(mul(m, tuple));
}

float3 rgbToXyz(float3 tuple) {
    const float3x3 m = float3x3(
        0.41239079926595948 , 0.35758433938387796, 0.18048078840183429 ,
        0.21263900587151036 , 0.71516867876775593, 0.072192315360733715,
        0.019330818715591851, 0.11919477979462599, 0.95053215224966058 
    );
    return mul(m, hsluv_toLinear(tuple));
}

float3 xyzToLuv(float3 tuple){
    float X = tuple.x;
    float Y = tuple.y;
    float Z = tuple.z;

    float L = hsluv_yToL(Y);
    
    float div = 1./dot(tuple,float3(1,15,3)); 

    return float3(
        1.,
        (52. * (X*div) - 2.57179),
        (117.* (Y*div) - 6.08816)
    ) * L;
}


float3 luvToXyz(float3 tuple) {
    float L = tuple.x;

    float U = tuple.y / (13.0 * L) + 0.19783000664283681;
    float V = tuple.z / (13.0 * L) + 0.468319994938791;

    float Y = hsluv_lToY(L);
    float X = 2.25 * U * Y / V;
    float Z = (3./V - 5.)*Y - (X/3.);

    return float3(X, Y, Z);
}

float3 luvToLch(float3 tuple) {
    float L = tuple.x;
    float U = tuple.y;
    float V = tuple.z;

    float C = length(tuple.yz);
    float H = degrees(atan2(V, U));
    if (H < 0.0) {
        H = 360.0 + H;
    }
    
    return float3(L, C, H);
}

float3 lchToLuv(float3 tuple) {
    float hrad = radians(tuple.b);
    return float3(
        tuple.r,
        cos(hrad) * tuple.g,
        sin(hrad) * tuple.g
    );
}

float3 hsluvToLch(float3 tuple) {
    tuple.g *= hsluv_maxChromaForLH(tuple.b, tuple.r) * .01;
    return tuple.bgr;
}

float3 lchToHsluv(float3 tuple) {
    tuple.g /= hsluv_maxChromaForLH(tuple.r, tuple.b) * .01;
    return tuple.bgr;
}

float3 hpluvToLch(float3 tuple) {
    tuple.g *= hsluv_maxSafeChromaForL(tuple.b) * .01;
    return tuple.bgr;
}

float3 lchToHpluv(float3 tuple) {
    tuple.g /= hsluv_maxSafeChromaForL(tuple.r) * .01;
    return tuple.bgr;
}

float3 lchToRgb(float3 tuple) {
    return xyzToRgb(luvToXyz(lchToLuv(tuple)));
}

float3 rgbToLch(float3 tuple) {
    return luvToLch(xyzToLuv(rgbToXyz(tuple)));
}

float3 hsluvToRgb(float3 tuple) {
    return lchToRgb(hsluvToLch(tuple));
}

float3 rgbToHsluv(float3 tuple) {
    return lchToHsluv(rgbToLch(tuple));
}

float3 hpluvToRgb(float3 tuple) {
    return lchToRgb(hpluvToLch(tuple));
}

float3 rgbToHpluv(float3 tuple) {
    return lchToHpluv(rgbToLch(tuple));
}

float3 luvToRgb(float3 tuple){
    return xyzToRgb(luvToXyz(tuple));
}

// allow float4's
float4   xyzToRgb(float4 c) {return float4(   xyzToRgb( float3(c.x,c.y,c.z) ), c.a);}
float4   rgbToXyz(float4 c) {return float4(   rgbToXyz( float3(c.x,c.y,c.z) ), c.a);}
float4   xyzToLuv(float4 c) {return float4(   xyzToLuv( float3(c.x,c.y,c.z) ), c.a);}
float4   luvToXyz(float4 c) {return float4(   luvToXyz( float3(c.x,c.y,c.z) ), c.a);}
float4   luvToLch(float4 c) {return float4(   luvToLch( float3(c.x,c.y,c.z) ), c.a);}
float4   lchToLuv(float4 c) {return float4(   lchToLuv( float3(c.x,c.y,c.z) ), c.a);}
float4 hsluvToLch(float4 c) {return float4( hsluvToLch( float3(c.x,c.y,c.z) ), c.a);}
float4 lchToHsluv(float4 c) {return float4( lchToHsluv( float3(c.x,c.y,c.z) ), c.a);}
float4 hpluvToLch(float4 c) {return float4( hpluvToLch( float3(c.x,c.y,c.z) ), c.a);}
float4 lchToHpluv(float4 c) {return float4( lchToHpluv( float3(c.x,c.y,c.z) ), c.a);}
float4   lchToRgb(float4 c) {return float4(   lchToRgb( float3(c.x,c.y,c.z) ), c.a);}
float4   rgbToLch(float4 c) {return float4(   rgbToLch( float3(c.x,c.y,c.z) ), c.a);}
float4 hsluvToRgb(float4 c) {return float4( hsluvToRgb( float3(c.x,c.y,c.z) ), c.a);}
float4 rgbToHsluv(float4 c) {return float4( rgbToHsluv( float3(c.x,c.y,c.z) ), c.a);}
float4 hpluvToRgb(float4 c) {return float4( hpluvToRgb( float3(c.x,c.y,c.z) ), c.a);}
float4 rgbToHpluv(float4 c) {return float4( rgbToHpluv( float3(c.x,c.y,c.z) ), c.a);}
float4   luvToRgb(float4 c) {return float4(   luvToRgb( float3(c.x,c.y,c.z) ), c.a);}
// allow 3 floats
float3   xyzToRgb(float x, float y, float z) {return   xyzToRgb( float3(x,y,z) );}
float3   rgbToXyz(float x, float y, float z) {return   rgbToXyz( float3(x,y,z) );}
float3   xyzToLuv(float x, float y, float z) {return   xyzToLuv( float3(x,y,z) );}
float3   luvToXyz(float x, float y, float z) {return   luvToXyz( float3(x,y,z) );}
float3   luvToLch(float x, float y, float z) {return   luvToLch( float3(x,y,z) );}
float3   lchToLuv(float x, float y, float z) {return   lchToLuv( float3(x,y,z) );}
float3 hsluvToLch(float x, float y, float z) {return hsluvToLch( float3(x,y,z) );}
float3 lchToHsluv(float x, float y, float z) {return lchToHsluv( float3(x,y,z) );}
float3 hpluvToLch(float x, float y, float z) {return hpluvToLch( float3(x,y,z) );}
float3 lchToHpluv(float x, float y, float z) {return lchToHpluv( float3(x,y,z) );}
float3   lchToRgb(float x, float y, float z) {return   lchToRgb( float3(x,y,z) );}
float3   rgbToLch(float x, float y, float z) {return   rgbToLch( float3(x,y,z) );}
float3 hsluvToRgb(float x, float y, float z) {return hsluvToRgb( float3(x,y,z) );}
float3 rgbToHsluv(float x, float y, float z) {return rgbToHsluv( float3(x,y,z) );}
float3 hpluvToRgb(float x, float y, float z) {return hpluvToRgb( float3(x,y,z) );}
float3 rgbToHpluv(float x, float y, float z) {return rgbToHpluv( float3(x,y,z) );}
float3   luvToRgb(float x, float y, float z) {return   luvToRgb( float3(x,y,z) );}
// allow 4 floats
float4   xyzToRgb(float x, float y, float z, float a) {return   xyzToRgb( float4(x,y,z,a) );}
float4   rgbToXyz(float x, float y, float z, float a) {return   rgbToXyz( float4(x,y,z,a) );}
float4   xyzToLuv(float x, float y, float z, float a) {return   xyzToLuv( float4(x,y,z,a) );}
float4   luvToXyz(float x, float y, float z, float a) {return   luvToXyz( float4(x,y,z,a) );}
float4   luvToLch(float x, float y, float z, float a) {return   luvToLch( float4(x,y,z,a) );}
float4   lchToLuv(float x, float y, float z, float a) {return   lchToLuv( float4(x,y,z,a) );}
float4 hsluvToLch(float x, float y, float z, float a) {return hsluvToLch( float4(x,y,z,a) );}
float4 lchToHsluv(float x, float y, float z, float a) {return lchToHsluv( float4(x,y,z,a) );}
float4 hpluvToLch(float x, float y, float z, float a) {return hpluvToLch( float4(x,y,z,a) );}
float4 lchToHpluv(float x, float y, float z, float a) {return lchToHpluv( float4(x,y,z,a) );}
float4   lchToRgb(float x, float y, float z, float a) {return   lchToRgb( float4(x,y,z,a) );}
float4   rgbToLch(float x, float y, float z, float a) {return   rgbToLch( float4(x,y,z,a) );}
float4 hsluvToRgb(float x, float y, float z, float a) {return hsluvToRgb( float4(x,y,z,a) );}
float4 rgbToHslul(float x, float y, float z, float a) {return rgbToHsluv( float4(x,y,z,a) );}
float4 hpluvToRgb(float x, float y, float z, float a) {return hpluvToRgb( float4(x,y,z,a) );}
float4 rgbToHpluv(float x, float y, float z, float a) {return rgbToHpluv( float4(x,y,z,a) );}
float4   luvToRgb(float x, float y, float z, float a) {return   luvToRgb( float4(x,y,z,a) );}

/*
END HSLUV-GLSL
*/

// Обертки для Shader Graph


// XYZ <-> RGB
void XYZToRGB_float(float3 xyz, out float3 rgb)
{
    rgb = xyzToRgb(xyz);
}
void RGBToXYZ_float(float3 rgb, out float3 xyz)
{
    xyz = rgbToXyz(rgb);
}

// XYZ <-> Luv
void XYZToLuv_float(float3 xyz, out float3 luv)
{
    luv = xyzToLuv(xyz);
}
void LuvToXYZ_float(float3 luv, out float3 xyz)
{
    xyz = luvToXyz(luv);
}

// Luv <-> LCH
void LuvToLCH_float(float3 luv, out float3 lch)
{
    lch = luvToLch(luv);
}
void LCHToLuv_float(float3 lch, out float3 luv)
{
    luv = lchToLuv(lch);
}

// HSLuv <-> LCH
void HSLuvToLCH_float(float3 hsluv, out float3 lch)
{
    lch = hsluvToLch(hsluv);
}
void LCHToHSLuv_float(float3 lch, out float3 hsluv)
{
    hsluv = lchToHsluv(lch);
}

// HPLuv <-> LCH
void HPLuvToLCH_float(float3 hpluv, out float3 lch)
{
    lch = hpluvToLch(hpluv);
}
void LCHToHPLuv_float(float3 lch, out float3 hpluv)
{
    hpluv = lchToHpluv(lch);
}

// LCH <-> RGB
void LCHToRGB_float(float3 lch, out float3 rgb)
{
    rgb = lchToRgb(lch);
}
void RGBToLCH_float(float3 rgb, out float3 lch)
{
    lch = rgbToLch(rgb);
}

// HSLuv/HPLuv <-> RGB
void HSLuvToRGB_float(float3 hsluv, out float3 rgb)
{
    rgb = hsluvToRgb(hsluv);
}
void HPLuvToRGB_float(float3 hpluv, out float3 rgb)
{
    rgb = hpluvToRgb(hpluv);
}
void RGBToHSLuv_float(float3 rgb, out float3 hsluv)
{
    hsluv = rgbToHsluv(rgb);
}
void RGBToHPLuv_float(float3 rgb, out float3 hpluv)
{
    hpluv = rgbToHpluv(rgb);
}


// XYZ <-> RGB
void XYZToRGB_float(float4 XYZ, out float4 RGB)
{
    RGB = xyzToRgb(XYZ);
}
void RGBToXYZ_float(float4 RGB, out float4 XYZ)
{
    XYZ = rgbToXyz(RGB);
}

// XYZ <-> Luv
void XYZToLuv_float(float4 XYZ, out float4 Luv)
{
    Luv = xyzToLuv(XYZ);
}
void LuvToXYZ_float(float4 Luv, out float4 XYZ)
{
    XYZ = luvToXyz(Luv);
}

// Luv <-> LCH
void LuvToLCH_float(float4 Luv, out float4 LCH)
{
    LCH = luvToLch(Luv);
}
void LCHToLuv_float(float4 LCH, out float4 Luv)
{
    Luv = lchToLuv(LCH);
}

// HSLuv <-> LCH
void HSLuvToLCH_float(float4 HSLuv, out float4 LCH)
{
    LCH = hsluvToLch(HSLuv);
}
void LCHToHSLuv_float(float4 LCH, out float4 HSLuv)
{
    HSLuv = lchToHsluv(LCH);
}

// HPLuv <-> LCH
void HPLuvToLCH_float(float4 HPLuv, out float4 LCH)
{
    LCH = hpluvToLch(HPLuv);
}
void LCHToHPLuv_float(float4 LCH, out float4 HPLuv)
{
    HPLuv = lchToHpluv(LCH);
}

// LCH <-> RGB
void LCHToRGB_float(float4 LCH, out float4 RGB)
{
    RGB = lchToRgb(LCH);
}
void RGBToLCH_float(float4 RGB, out float4 LCH)
{
    LCH = rgbToLch(RGB);
}

// HSLuv <-> RGB
void HSLuvToRGB_float(float4 HSLuv, out float4 RGB)
{
    RGB = hsluvToRgb(HSLuv);
}
void RGBToHSLuv_float(float4 RGB, out float4 HSLuv)
{
    HSLuv = rgbToHsluv(RGB);
}

// HPLuv <-> RGB
void HPLuvToRGB_float(float4 HPLuv, out float4 RGB)
{
    RGB = hpluvToRgb(HPLuv);
}
void RGBToHPLuv_float(float4 RGB, out float4 HPLuv)
{
    HPLuv = rgbToHpluv(RGB);
}

// Luv <-> RGB
void LuvToRGB_float(float4 Luv, out float4 RGB)
{
    RGB = luvToRgb(Luv);
}

// XYZ <-> RGB
void XYZToRGB_Components(float x, float y, float z, float a, out float4 RGB)
{
    RGB = xyzToRgb(float4(x, y, z, a));
}
void RGBToXYZ_Components(float r, float g, float b, float a, out float4 XYZ)
{
    XYZ = rgbToXyz(float4(r, g, b, a));
}

// XYZ <-> Luv
void XYZToLuv_Components(float x, float y, float z, float a, out float4 Luv)
{
    Luv = xyzToLuv(float4(x, y, z, a));
}
void LuvToXYZ_Components(float u, float v, float w, float a, out float4 XYZ)
{
    XYZ = luvToXyz(float4(u, v, w, a));
}

// Luv <-> LCH
void LuvToLCH_Components(float u, float v, float w, float a, out float4 LCH)
{
    LCH = luvToLch(float4(u, v, w, a));
}
void LCHToLuv_Components(float l, float c, float h, float a, out float4 Luv)
{
    Luv = lchToLuv(float4(l, c, h, a));
}

// HSLuv <-> LCH
void HSLuvToLCH_Components(float h, float s, float l, float a, out float4 LCH)
{
    LCH = hsluvToLch(float4(h, s, l, a));
}
void LCHToHSLuv_Components(float l, float c, float h, float a, out float4 HSLuv)
{
    HSLuv = lchToHsluv(float4(l, c, h, a));
}

// HPLuv <-> LCH
void HPLuvToLCH_Components(float h, float p, float l, float a, out float4 LCH)
{
    LCH = hpluvToLch(float4(h, p, l, a));
}
void LCHToHPLuv_Components(float l, float c, float h, float a, out float4 HPLuv)
{
    HPLuv = lchToHpluv(float4(l, c, h, a));
}

// LCH <-> RGB
void LCHToRGB_Components(float l, float c, float h, float a, out float4 RGB)
{
    RGB = lchToRgb(float4(l, c, h, a));
}
void RGBToLCH_Components(float r, float g, float b, float a, out float4 LCH)
{
    LCH = rgbToLch(float4(r, g, b, a));
}

// HSLuv <-> RGB
void HSLuvToRGB_Components(float h, float s, float l, float a, out float4 RGB)
{
    RGB = hsluvToRgb(float4(h, s, l, a));
}
void RGBToHSLuv_Components(float r, float g, float b, float a, out float4 HSLuv)
{
    HSLuv = rgbToHsluv(float4(r, g, b, a));
}

// HPLuv <-> RGB
void HPLuvToRGB_Components(float h, float p, float l, float a, out float4 RGB)
{
    RGB = hpluvToRgb(float4(h, p, l, a));
}
void RGBToHPLuv_Components(float r, float g, float b, float a, out float4 HPLuv)
{
    HPLuv = rgbToHpluv(float4(r, g, b, a));
}

// Luv <-> RGB
void LuvToRGB_Components(float l, float u, float v, float a, out float4 RGB)
{
    RGB = luvToRgb(float4(l, u, v, a));
}


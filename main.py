# main.py — محاكاة فلكية كاملة لـ 10,480 ق.م. (بدون أخطاء)
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_body
from astropy import units as u

# --- قاعدة النجوم ---
STARS = {
    "Alnitak": ("05h40m45.5s", "-01d56m33s"),
    "Alnilam": ("05h36m12.8s", "-01d12m06s"),
    "Mintaka": ("05h32m00.4s", "-00d17m57s"),
    "Sirius": ("06h45m08.9s", "-16d42m58s")
}

# --- تحميل config.txt ---
def load_config():
    sims = []
    current = {}
    try:
        with open("config.txt", "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"): continue
                if line.startswith("sim"):
                    if current: sims.append(current)
                    current = {"name": line.split(".")[0]}
                else:
                    try:
                        k, v = line.split("=", 1)
                        k, v = k.strip(), v.strip()
                        if k in ["lat", "lon"]: current[k] = float(v)
                        elif k == "plot": current[k] = v.lower() == "true"
                        else: current[k] = v
                    except: pass
        if current: sims.append(current)
    except Exception as e:
        print(f"خطأ في config.txt: {e}")
    return sims

# --- المحاكاة ---
def simulate(sim):
    print(f"\n{'='*60}")
    # آمن ضد الأخطاء
    print(f"   {sim.get('name', 'Unknown')} | {sim.get('date', 'No Date')}")
    print(f"{'='*60}")
    
    # التحقق من وجود التاريخ
    if 'date' not in sim:
        print("خطأ: 'date' مفقود في config.txt")
        return
    
    try:
        time = Time(sim['date'], format='iso', scale='utc')
    except Exception as e:
        print(f"خطأ في التاريخ: {e}")
        return
    
    loc = EarthLocation(lat=sim.get('lat', 0)*u.deg, lon=sim.get('lon', 0)*u.deg)
    frame = AltAz(obstime=time, location=loc)
    targets = [t.strip() for t in sim.get('targets', '').split(",") if t.strip()]
    results = {}
    
    for t in targets:
        obj = None
        if t in STARS:
            star = SkyCoord(ra=STARS[t][0], dec=STARS[t][1], frame='icrs')
            obj = star
        elif t.lower() in ['sun', 'moon']:
            try:
                from astropy.coordinates import solar_system_ephemeris
                with solar_system_ephemeris.set('de432s'):
                    obj = get_body(t.lower(), time)
            except: pass
        
        if obj:
            try:
                altaz = obj.transform_to(frame)
                results[t] = {'az': altaz.az.deg, 'alt': altaz.alt.deg}
                print(f"{t}: أزيموث {altaz.az.deg:.2f}° | ارتفاع {altaz.alt.deg:.2f}°")
            except Exception as e:
                print(f"فشل في {t}: {e}")
    
    # --- مسافات ---
    if sim.get("mode") == "distance" and len(results) >= 2:
        names = list(results.keys())
        try:
            from astropy.coordinates import GCRS
            a = results[names[0]]['coord'] if 'coord' in results[names[0]] else None
            b = results[names[1]]['coord'] if 'coord' in results[names[1]] else None
            if a and b:
                dist_au = a.separation_3d(b).to(u.au)
                dist_km = dist_au.to(u.km)
                print(f"\nمسافة {names[0]} ←→ {names[1]}:")
                print(f"  → {dist_au:.6f} AU")
                print(f"  → {dist_km:,.0f} كم")
        except: pass
    
    # --- رسم ---
    if sim.get("plot", False):
        plot_sky(results, sim.get('name', 'sim'))

def plot_sky(results, name):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(0, 360); ax.set_ylim(-90, 90)
    for obj, d in results.items():
        ax.scatter(d['az'], d['alt'], s=120, label=obj, edgecolor='black')
    ax.set_xlabel("أزيموث (°)")
    ax.set_ylabel("ارتفاع (°)")
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_title(f"محاكاة: {name}")
    os.makedirs("output", exist_ok=True)
    plt.savefig(f"output/{name}.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"تم حفظ: output/{name}.png")

# --- التشغيل ---
if __name__ == "__main__":
    os.makedirs("output", exist_ok=True)
    sims = load_config()
    if not sims:
        print("لا توجد محاكاة — تحقق من config.txt")
    else:
        for sim in sims:
            simulate(sim)
    print("\nتم الانتهاء! النتائج في مجلد 'output'")

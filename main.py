

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_moon
from astropy.time import Time
import astropy.units as u

# =========================================================
# 🔹 تحميل الإعدادات من config.txt
# =========================================================
def load_config():
    sims = {}
    with open("config.txt", "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("---"):
                continue
            if "=" not in line:
                continue
            left, right = line.split("=", 1)
            left = left.strip()
            val = right.strip()
            if "." not in left:
                continue
            sim_id, key = left.split(".", 1)
            if sim_id not in sims:
                sims[sim_id] = {}
            if key in ("lat", "lon"):
                sims[sim_id][key] = float(val)
            elif key == "plot":
                sims[sim_id][key] = val.lower() in ("true", "1", "yes", "y")
            else:
                sims[sim_id][key] = val

    def order_key(k):
        m = re.search(r"(\d+)$", k)
        return int(m.group(1)) if m else 0

    return [sims[k] for k in sorted(sims.keys(), key=order_key)]

# =========================================================
# 🔹 التحقق من صحة الإدخال
# =========================================================
def validate_sim(sim):
    required_keys = ["name", "lat", "lon", "date", "targets", "plot", "mode"]
    for key in required_keys:
        if key not in sim:
            print(f"❌ Missing key '{key}' in simulation: {sim.get('name', 'Unnamed')}")
            return False
    try:
        float(sim["lat"])
        float(sim["lon"])
        parse_time(sim["date"])
    except Exception as e:
        print(f"❌ Invalid value in simulation '{sim['name']}': {e}")
        return False
    if sim["mode"].lower() not in ("alignment", "distance", "full"):
        print(f"❌ Unknown mode '{sim['mode']}' in simulation '{sim['name']}'")
        return False
    return True

# =========================================================
# 🔹 دالة لتحليل التاريخ (قبل الميلاد)
# =========================================================
def parse_time(date_str):
    if not date_str:
        raise ValueError("Date missing in config.txt")
    try:
        return Time(date_str, format="iso", scale="utc")
    except Exception:
        pass
    m = re.match(r"^(-?\d+)", date_str)
    if m:
        return Time(int(m.group(1)), format="byear")
    raise ValueError(f"Unrecognized date format: {date_str}")

# =========================================================
# 🔹 محاكاة فلكية
# =========================================================
def simulate(sim):
    name = sim.get("name", "Unnamed")
    mode = sim.get("mode", "alignment").lower()
    print(f"\n=== Running {name} [{mode}] ===")

    try:
        t = parse_time(sim.get("date", ""))
    except Exception as e:
        print(f"❌ Error parsing date: {e}")
        return

    try:
        loc = EarthLocation(lat=sim["lat"] * u.deg, lon=sim["lon"] * u.deg)
    except Exception as e:
        print(f"❌ Invalid location: {e}")
        return

    altaz_frame = AltAz(obstime=t, location=loc)
    targets = [t.strip().lower() for t in sim["targets"].split(",")]
    coords = {}

    for target in targets:
        try:
            if target in ["sun", "☉"]:
                obj = get_sun(t).transform_to(altaz_frame)
            elif target in ["moon", "🌙"]:
                obj = get_moon(t).transform_to(altaz_frame)
            elif target in ["alnitak", "alnilam", "mintaka"]:
                stars = {
                    "alnitak": SkyCoord.from_name("Alnitak"),
                    "alnilam": SkyCoord.from_name("Alnilam"),
                    "mintaka": SkyCoord.from_name("Mintaka"),
                }
                obj = stars[target].transform_to(altaz_frame)
            elif target in ["sirius", "α canis majoris"]:
                obj = SkyCoord.from_name("Sirius").transform_to(altaz_frame)
            else:
                print(f"⚠️ Unknown target: {target}")
                continue
            coords[target] = (obj.alt.degree, obj.az.degree)
        except Exception as e:
            print(f"❌ Error with target '{target}': {e}")

    if mode == "distance":
        for target in targets:
            if target == "sun":
                dist = get_sun(t).distance.to(u.km)
                print(f"   ☉ Sun Distance = {dist:.2f}")
            elif target == "moon":
                dist = get_moon(t).distance.to(u.km)
                print(f"   🌙 Moon Distance = {dist:.2f}")
        return

    if sim.get("plot", False) and coords:
        if not os.path.exists("output"):
            os.makedirs("output")
        plt.figure(figsize=(8, 5))
        for k, (alt, az) in coords.items():
            plt.scatter(az, alt, label=k.capitalize())
        plt.xlabel("Azimuth (°)")
        plt.ylabel("Altitude (°)")
        plt.title(f"{name} — {t.iso}")
        plt.legend()
        plt.grid(True)
        outfile = os.path.join("output", f"{name}.png")
        plt.savefig(outfile, dpi=150)
        plt.close()
        print(f"✅ Saved plot → {outfile}")

    for k, (alt, az) in coords.items():
        print(f"   {k.capitalize():<10} Alt={alt:>7.2f}°,  Az={az:>7.2f}°")

    if mode in ("alignment", "full") and len(coords) >= 2:
        azimuths = [az for _, az in coords.values()]
        spread = max(azimuths) - min(azimuths)
        print(f"   🔸 Azimuth spread = {spread:.2f}°")

# =========================================================
# 🔹 تشغيل البرنامج
# =========================================================
if __name__ == "__main__":
    sims = load_config()
    for s in sims:
        if validate_sim(s):
            simulate(s)
        else:
            print("⚠️ Skipping invalid simulation.\n")
    print("\n✅ Simulation Completed.")

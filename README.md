# WHITE TIGER · Projectile Calculator

เครื่องมือจำลองวิถีกระสุนแบบ Interactive สำหรับคำนวณมุมยิง, เวลาหน่วง (Delay Time) และค่าทางไฟฟ้าที่เกี่ยวข้อง พร้อมแสดงผลกราฟวิถีกระสุนแบบ Real-time บน Canvas 2D

---

## รายละเอียดโปรเจกต์ (Description)

**WHITE TIGER Projectile Calculator** เป็นเว็บแอปพลิเคชันที่ออกแบบมาเพื่อช่วยในการคำนวณพารามิเตอร์ของเครื่องยิงกระสุน (Projectile Launcher) โดยมีความสามารถหลัก ดังนี้

- **คำนวณมุมยิง (Launch Angle)** — ใช้สมการเชิงประจักษ์ (Empirical Quadratic) เพื่อหามุมที่เหมาะสมจากระยะเป้าหมาย
- **คำนวณเวลาหน่วง (Delay Time)** — คำนวณจากมุมเซ็นเซอร์, ความเร็วรอบ (RPM) และจำนวนรอบการหมุน
- **คำนวณค่าทางไฟฟ้า** — หาค่าความต้านทาน (R) และแรงดันไฟฟ้า (V) ที่ต้องการ
- **กราฟวิถีกระสุน (Trajectory Graph)** — แสดงวิถีกระสุนแบบ Interactive บน Canvas 2D พร้อมตำแหน่งเป้าหมาย
- **ระบบ Preset History** — บันทึกและโหลดค่าที่เคยคำนวณไว้ผ่าน LocalStorage

---

## โครงสร้างไฟล์ (Project Structure)

```
M3.1/
├── index.html    # โครงสร้าง HTML หลัก (UI Layout 3 คอลัมน์)
├── styles.css    # ระบบออกแบบ CSS (CSS Variables, Responsive, Animation)
├── main.js       # ลอจิกทั้งหมด (Physics, Timing, Canvas, UI Controller)
└── README.md     # เอกสารอธิบายโปรเจกต์
```

---

## ฟีเจอร์หลัก (Features)

| ฟีเจอร์ | รายละเอียด |
| --- | --- |
| **Target Calculator** | กรอกระยะวัดเป้า (y) แล้วระยะเป้าหมาย (sₓ) จะคำนวณอัตโนมัติ |
| **Launch Angle Solver** | แก้สมการกำลังสองหามุมยิง (High Arc) ในช่วง 46°–73° |
| **Trajectory Graph** | แสดงวิถีกระสุนบน Canvas 2D พร้อมเส้นอ้างอิงความสูง |
| **Angle Adjuster** | ปรับมุมด้วยปุ่ม +/− แล้วกราฟอัปเดตแบบ Real-time |
| **Delay Time Calculator** | คำนวณเวลาหน่วงจากมุมเซ็นเซอร์, RPM และจำนวนรอบ |
| **Electrical Data** | คำนวณค่า Resistance (kΩ) และ Voltage (V) |
| **Dirty Tracking** | ระบบตรวจจับการเปลี่ยนแปลงค่า พร้อมแจ้งเตือนให้กด CALCULATE |
| **Preset History** | บันทึกค่าอัตโนมัติ สูงสุด 10 รายการ เรียกใช้ซ้ำได้ |
| **Responsive Design** | รองรับหน้าจอทุกขนาด ตั้งแต่ Desktop ถึง Mobile |

---

## ค่าคงที่ในระบบ (Design Constants)

| พารามิเตอร์ | ค่า | หมายเหตุ |
| --- | --- | --- |
| แรงโน้มถ่วง (g) | 9.81 m/s² | ค่ามาตรฐาน |
| ความสูงเครื่องยิง (y₀) | 380 mm | จุดปล่อยกระสุน |
| ความสูงเป้าหมาย (yₜ) | 430 mm | ตำแหน่งเป้า |
| ช่วงมุมที่ถูกต้อง (θ) | 46°–73° | ขอบเขตการคำนวณ |

---

## สถาปัตยกรรมระบบ (Architecture)

โปรเจกต์ใช้สถาปัตยกรรมแบบ **OOP (Object-Oriented Programming)** โดยแบ่งเป็นโมดูลอิสระ:

```
User Input
  → InputController (จัดการ Input ทั้งหมด)
  → PhysicsEngine (คำนวณมุมยิง)
  → AngleController (จัดการสถานะมุม)
  → CanvasRenderer (วาดกราฟ)
  → UIRenderer (อัปเดต DOM)
  → TimingController → TimingEngine (คำนวณเวลาหน่วง)
  → HistoryController (จัดการประวัติ)
```
---

## หมายเหตุ (Notes)

- สมการที่ใช้ในการคำนวณเป็นสมการเชิงประจักษ์ (Empirical) ที่ได้จากการทดสอบจริง
- ข้อมูลประวัติจัดเก็บใน `localStorage` ของเบราว์เซอร์ ไม่มีการส่งข้อมูลไปยัง Server
- โปรเจกต์นี้เป็น **Static Web Application** ไม่ต้องใช้ Backend หรือ Database

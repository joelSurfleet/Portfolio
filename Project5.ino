#include <Wire.h>
#include <Adafruit_Sensor.h>
#include <Adafruit_BNO055.h>
#include <utility/imumaths.h>
  
Adafruit_BNO055 bno = Adafruit_BNO055(4);

const int LEDred1 = 2;
const int LEDred2 = 3;
const int LEDred3 = 4;
const int LEDgre1 = 6;
const int LEDgre2 = 7;
const int LEDgre3 = 8;
const int L = 10;
const int R = 11;

float time = 250;

void setup() {
  pinMode(LEDred1,OUTPUT);
  pinMode(LEDred2,OUTPUT);
  pinMode(LEDred3,OUTPUT);
  pinMode(LEDgre1,OUTPUT);
  pinMode(LEDgre2,OUTPUT);
  pinMode(LEDgre3,OUTPUT);

  Serial.begin(9600);
  Serial.println("Orientation Sensor Test"); Serial.println("");
  
  if(!bno.begin()) {
    Serial.print("Ooops, no BNO055 detected ... Check your wiring or I2C ADDR!");
    while(1);
  }
  
  delay(1000);
    
  bno.setExtCrystalUse(true);
  digitalWrite(L,HIGH);
  digitalWrite(R,HIGH);
}

void loop() {
  sensors_event_t event; 
  bno.getEvent(&event);

  // Display the floating point data 
  Serial.print("X: ");
  Serial.print(event.orientation.x, 4);
  Serial.print("\tY: ");
  Serial.print(event.orientation.y, 4);
  Serial.print("\tZ: ");
  Serial.print(event.orientation.z, 4);
  Serial.println("");
  
  if(event.orientation.z < 50 && event.orientation.z > 0){
    move((abs(event.orientation.z)/50)*100,1);
  }
  if(event.orientation.z > -50 && event.orientation.z < 0){
    move((abs(event.orientation.z)/50)*100,0);
  }
  else {
    move(100,2);
  }
  if(event.orientation.z > 30) {
    digitalWrite(LEDgre1,HIGH);
  }
  else {
    digitalWrite(LEDgre1,LOW);
  }
  if(event.orientation.z > 40) {
    digitalWrite(LEDgre2,HIGH);
  }
  else {
    digitalWrite(LEDgre2,LOW);
  }
  if(event.orientation.z > 50) {
    digitalWrite(LEDgre3,HIGH);
  }
  else {
    digitalWrite(LEDgre3,LOW);
  }
  if(event.orientation.z < -30) {
    digitalWrite(LEDred1,HIGH);
  }
  else {
    digitalWrite(LEDred1,LOW);
  }
  if(event.orientation.z < -40) {
    digitalWrite(LEDred2,HIGH);
  }
  else {
    digitalWrite(LEDred2,LOW);
  }
  if(event.orientation.z < -50) {
    digitalWrite(LEDred3,HIGH);
  }
  else {
    digitalWrite(LEDred3,LOW);
  }
  delay(10);
}

void move(float speed, int _direction)
{
// inputes speed as a percent
// needs global variable "time" which is the time constant of the motor

switch(_direction){
  case 1:
    digitalWrite(L,HIGH);
    digitalWrite(R,HIGH);
    delayMicroseconds(speed*time);
    digitalWrite(L,LOW);
    delayMicroseconds((100-speed)*time);
    break;
  case 0:
    digitalWrite(L,HIGH);
    digitalWrite(R,HIGH);
    delayMicroseconds(speed*time);
    digitalWrite(R,LOW);
    delayMicroseconds((100-speed)*time);
    break;
  default:
    digitalWrite(L,HIGH);
    digitalWrite(R,HIGH);
    break;
  }
}

void stop() {
// some code to stop the motor by disabling the L293D
digitalWrite(L,0);
digitalWrite(R,0);
}
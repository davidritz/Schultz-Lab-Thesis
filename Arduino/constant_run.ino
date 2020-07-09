#define inPin0 0
#define inPin1 1
#define inPin2 2
#define outPin8 8
#define outPin11 11
#define outPin12 12

int pd0 = 0;
int pd1 = 1;
int pd2 = 2;
int pd8 = 8;
int pd11 = 11;
int pd12 = 12;

void setup() { 
  // set photodiodes to input. Will check voltage at A0, A1, A2
  pinMode(pd0, INPUT);
  pinMode(pd1, INPUT);
  pinMode(pd2, INPUT);
  pinMode(pd8, OUTPUT);
  pinMode(pd11, OUTPUT);
  pinMode(pd12, OUTPUT);
  Serial.begin(9600);
  Serial.println();
}
void loop() {
  // put your main code here, to run repeatedly:
  digitalWrite(8, LOW);
  digitalWrite(11, LOW);
  digitalWrite(12, LOW);
  digitalWrite(12, HIGH); // sets the digital pin 8 on
  delay(2000);            // waits for 2 seconds

  // assign pins
  int pinRead0 = analogRead(inPin0);
  int pinRead1 = analogRead(inPin1);
  int pinRead2 = analogRead(inPin2);
  
  // calculate V
  float pVolt0 = pinRead0;
  float pVolt1 = pinRead1;
  float pVolt2 = pinRead2;
  
  float GFP = pVolt0;
  float Laser = pVolt1;
  float mCherry = pVolt2;
  // print in terminal
  Serial.print(" -GFP-> ");
  Serial.println(GFP);

  Serial.print(" -Laser-> ");
  Serial.println(Laser);

  Serial.print(" -mCherry-> ");
  Serial.println(mCherry);
  Serial.print(" ");
  Serial.println("----------");
}

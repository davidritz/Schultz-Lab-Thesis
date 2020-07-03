// For this code, I used my phone's flashlight on each photodiode to make sure voltage would increase correctly

// assign and call photodiode pins

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
  digitalWrite(8, LOW);
  digitalWrite(11, LOW);
  digitalWrite(12, LOW);
  digitalWrite(8, HIGH); // sets the digital pin 8 on
  delay(10000);            // waits for 2 seconds

  // assign pins
  int pinRead0 = analogRead(inPin0);
  int pinRead1 = analogRead(inPin1);
  int pinRead2 = analogRead(inPin2);
  
  // calculate V
  float pVolt0 = pinRead0 / 80.0 * 5.0;
  float pVolt1 = pinRead1 / 80.0 * 5.0;
  float pVolt2 = pinRead2 / 80.0 * 5.0;
  
  float GFP = pVolt0;
  float Laser = pVolt1;
  float mCherry = pVolt2;
  
  for (int i=1; i<10; i++)
  {
  pVolt0 = pinRead0 / 80.0 * 5.0;
  GFP = GFP+pVolt0;
  pVolt1 = pinRead1 / 80.0 * 5.0;
  Laser = Laser+pVolt1;
  pVolt2 = pinRead2 / 80.0 * 5.0;
  mCherry = mCherry+pVolt2;
  delay(2000);
  }
  
  GFP = GFP/10;
  Laser = Laser/10;
  mCherry = mCherry/10;
  
  // print in terminal
  Serial.print(" -GFP-> ");
  Serial.println(GFP);

  Serial.print(" -Laser-> ");
  Serial.println(Laser);

  Serial.print(" -mCherry-> ");
  Serial.println(mCherry);
  Serial.print(" ");
  Serial.println("----------");

  digitalWrite(8, LOW);  // sets the digital pin 8 off
  delay(2000);            // waits for 2 seconds
  digitalWrite(11, HIGH);  // sets the digital pin 11 on
  delay(2000);            // waits for 2 seconds

  pinRead0 = analogRead(inPin0);
  pinRead1 = analogRead(inPin1);
  pinRead2 = analogRead(inPin2);
  
  // calculate V
  pVolt0 = pinRead0 / 80.0 * 5.0;
  pVolt1 = pinRead1 / 80.0 * 5.0;
  pVolt2 = pinRead2 / 80.0 * 5.0;
  
  GFP = pVolt0;
  Laser = pVolt1;
  mCherry = pVolt2;

  for (int i=1; i<10; i++)
  {
  pVolt0 = pinRead0 / 80.0 * 5.0;
  GFP = GFP+pVolt0;
  pVolt1 = pinRead1 / 80.0 * 5.0;
  Laser = Laser+pVolt1;
  pVolt2 = pinRead2 / 80.0 * 5.0;
  mCherry = mCherry+pVolt2;
  delay(2000);
  }
  
  GFP = GFP/10;
  Laser = Laser/10;
  mCherry = mCherry/10;
  
  // print in terminal
  Serial.print(" -GFP-> ");
  Serial.println(GFP);

  Serial.print(" -Laser-> ");
  Serial.println(Laser);

  Serial.print(" -mCherry-> ");
  Serial.println(mCherry);
  Serial.print(" ");
  Serial.println("----------");

  digitalWrite(11, LOW);  // sets the digital pin 11 off
  delay(2000);            // waits for 2 seconds
  digitalWrite(12, HIGH);  // sets the digital pin 12 on
  delay(2000);

  pinRead0 = analogRead(inPin0);
  pinRead1 = analogRead(inPin1);
  pinRead2 = analogRead(inPin2);
  
  // calculate V
  pVolt0 = pinRead0 / 80.0 * 5.0;
  pVolt1 = pinRead1 / 80.0 * 5.0;
  pVolt2 = pinRead2 / 80.0 * 5.0;
  
  GFP = pVolt0;
  Laser = pVolt1;
  mCherry = pVolt2;

  for (int i=1; i<10; i++)
  {
  pVolt0 = pinRead0 / 80.0 * 5.0;
  GFP = GFP+pVolt0;
  pVolt1 = pinRead1 / 80.0 * 5.0;
  Laser = Laser+pVolt1;
  pVolt2 = pinRead2 / 80.0 * 5.0;
  mCherry = mCherry+pVolt2;
  delay(2000);
  }
  
  GFP = GFP/10;
  Laser = Laser/10;
  mCherry = mCherry/10;
  
  // print in terminal
  Serial.print(" -GFP-> ");
  Serial.println(GFP);

  Serial.print(" -Laser-> ");
  Serial.println(Laser);

  Serial.print(" -mCherry-> ");
  Serial.println(mCherry);
  Serial.print(" ");
  Serial.println("----------");

  digitalWrite(12, LOW);  // sets the digital pin 12 off
  delay(2000);            // waits for 2 seconds
  
  delay(100);
}

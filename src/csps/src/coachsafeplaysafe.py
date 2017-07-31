from kivy.uix.screenmanager import ScreenManager, Screen, SwapTransition
from kivy.uix.listview import ListItemButton
from kivy.properties import ObjectProperty
from kivy.uix.gridlayout import GridLayout
from kivy.uix.pagelayout import PageLayout
from kivy.uix.boxlayout import BoxLayout
from kivy.properties import ListProperty
from kivy.app import App
import threespace_api as ts_api
from random import gauss
import sqlite3
import os


mnnno_sensor_api = False
no_sensor = False
notEstablished = True

if os.name == 'posix':
    no_sensor_api = True


class StudentListButton(ListItemButton):
    pass

class ScreenOne(Screen):
    title = 'Zero'
   
    fn = ObjectProperty()
    ln = ObjectProperty()
    student_list = ObjectProperty()

    def submit_student(self):
        student_name = self.fn.text +" "+ self.ln.text
        first_name = self.fn.text
        last_name =self.ln.text
        
        self.student_list.adapter.data.extend([student_name])
        self.student_list._trigger_reset_populate()
        
        con = sqlite3.connect('players.db')
        cur = con.cursor() 
        cur.execute("Create TABLE if not exists Players(first_name TEXT, last_name TEXT)")
        cur.execute("INSERT INTO Players (first_name, last_name) VALUES (?,?)", (first_name, last_name))
        con.commit()   
        cur.execute("SELECT * FROM Players")
        data = cur.fetchall()

    def delete_student(self, *args):
        if self.student_list.adapter.selection:
            selection = self.student_list.adapter.selection[0].text
            self.student_list.adapter.data.remove(selection)
            self.student_list._trigger_reset_populate()

        first_name = self.fn.text
        last_name =self.ln.text
        con = sqlite3.connect('players.db')
        cur = con.cursor() 
        cur.execute("SELECT * FROM Players")
        cur.execute('DELETE FROM Players WHERE first_name=(?) and last_name=(?)', (first_name, last_name))
        con.commit()


    def replace_student(self, *args):
        if self.student_list.adapter.selection:

            selection = self.student_list.adapter.selection[0].text
            self.student_list.adapter.data.remove(selection)
            student_name = self.fn.text +" "+ self.ln.text
            self.student_list.adapter.data.extend([student_name])
            self.student_list._trigger_reset_populate()
            print (selection)

        first_name = self.fn.text
        last_name =self.ln.text
        con = sqlite3.connect('players.db')
        cur = con.cursor() 
        cur.execute("SELECT * FROM Players")
        cur.execute('UPDATE Players SET first_name=(?) WHERE first_name=(?)', (first_name))
        con.commit()


class ScreenTwo(Screen):
    pass
class CoachSafePlaySafeApp(App):
    def build(self):
        screen_manager = ScreenManager(transition=SwapTransition())
        screen_manager.add_widget(ScreenOne(name="MPV"))
        screen_manager.add_widget(ScreenTwo(name="VidSensor"))
        return screen_manager


if __name__ == "__main__":
    CoachSafePlaySafeApp().run()
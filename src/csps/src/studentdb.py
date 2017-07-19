from kivy.app import App
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.pagelayout import PageLayout
import sqlite3
from kivy.uix.label import Label
from kivy.uix.textinput import TextInput
import random

from kivy.properties import ObjectProperty
from kivy.uix.listview import ListItemButton
import multiprocessing

class StudentListButton(ListItemButton):
    pass


class StudentDB(BoxLayout):

    fn = ObjectProperty()
    ln = ObjectProperty()
    student_list = ObjectProperty()


    def submit_student(self):
  
        # Get the student name from the TextInputs
        student_name = self.fn.text +" "+ self.ln.text
        first_name = self.fn.text
        last_name =self.ln.text

        # Add the student to the ListView
        self.student_list.adapter.data.extend([student_name])

        # Reset the ListView
        self.student_list._trigger_reset_populate()

        con = sqlite3.connect('playerscom.db')
        cur = con.cursor() 
        cur.execute("Create TABLE if not exists Players(first_name TEXT primary key, last_name TEXT unique)")
        cur.execute("INSERT INTO Players (first_name, last_name) VALUES (?,?)", (first_name, last_name))
        con.commit()   
        cur.execute("SELECT * FROM Players")
        data = cur.fetchall()

        # obj = CamApp()
        # p = multiprocessing.Process(target=obj.run)
        # p.start()

    def delete_student(self, *args):

        # If a list item is selected
        if self.student_list.adapter.selection:

            # Get the text from the item selected
            selection = self.student_list.adapter.selection[0].text

            # Remove the matching item
            self.student_list.adapter.data.remove(selection)

            # Reset the ListView
            self.student_list._trigger_reset_populate()

        first_name = self.fn.text
        last_name =self.ln.text

        con = sqlite3.connect('test.db')
        cur = con.cursor() 
        cur.execute("DELETE FROM Players WHERE first_name=self.fn.text")
        con.commit()   

 
    def replace_student(self, *args):

        # If a list item is selected
        if self.student_list.adapter.selection:

            # Get the text from the item selected
            selection = self.student_list.adapter.selection[0].text

            # Remove the matching item
            self.student_list.adapter.data.remove(selection)

            # Get the student name from the TextInputs
            student_name = self.fn.text +" "+ self.ln.text

            # Add the updated data to the list
            self.student_list.adapter.data.extend([student_name])

            # Reset the ListView
            self.student_list._trigger_reset_populate()

class StudentDBApp(App):
   
    def build(self):
        return StudentDB()

if __name__ == '__main__':
    StudentDBApp().run()

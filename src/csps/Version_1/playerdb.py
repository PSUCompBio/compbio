from kivy.app import App
from kivy.uix.boxlayout import BoxLayout
from kivy.properties import ObjectProperty
from kivy.uix.listview import ListItemButton
from tkinter import *
# from tkinter import tkk
# test git commit
import sqlite3

# class playerDB
#     db_conn =0  
#     cursor = 0
#     curr_player = 0

 
class PlayerListButton(ListItemButton):
    pass
 
class PlayerDB(BoxLayout):
 
    # Connects the value in the TextInput widget to these
    # fields
    player_name_text_input = ObjectProperty()
    team_text_input = ObjectProperty()
    player_list = ObjectProperty()
 
    def submit_player(self):
 
        # Get the student name from the TextInputs
        player_name = self.player_name_text_input.text + " " + self.team_text_input.text
 
        # Add the student to the ListView
        self.player_list.adapter.data.extend([player_name])
 
        # Reset the ListView
        self.player_list._trigger_reset_populate()
 
    def delete_player(self, *args):
 
        # If a list item is selected
        if self.player_list.adapter.selection:
 
            # Get the text from the item selected
            selection = self.player_list.adapter.selection[0].text
 
            # Remove the matching item
            self.player_list.adapter.data.remove(selection)
 
            # Reset the ListView
            self.player_list._trigger_reset_populate()
 
    def replace_player(self, *args):
 
        # If a list item is selected
        if self.player_list.adapter.selection:
 
            # Get the text from the item selected
            selection = self.player_list.adapter.selection[0].text
 
            # Remove the matching item
            self.player_list.adapter.data.remove(selection)
 
            # Get the student name from the TextInputs
            player_name = self.player_name_text_input.text + " " + self.team_text_input.text
 
            # Add the updated data to the list
            self.player_list.adapter.data.extend([player_name])
 
            # Reset the ListView
            self.player_list._trigger_reset_populate()

 
 
class PlayerDBApp(App):
    def build(self):
        return PlayerDB()
 
 
dbApp = PlayerDBApp()
 
dbApp.run()
 

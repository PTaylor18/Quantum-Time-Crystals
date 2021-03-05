using TerminalMenus

options = ["apple", "orange", "grape", "strawberry",
                "blueberry", "peach", "lemon", "lime", "Bananna",]

                # `pagesize` is the number of items to be displayed at a time.
                #  The UI will scroll if the number of options is greater
                #   than the `pagesize`
                menu = RadioMenu(options, pagesize=10)

                # `request` displays the menu and returns the index after the
                #   user has selected a choice
                choice = request("Choose your favorite fruit:", menu)

                if choice != -1
                    println("Your favorite fruit is ", options[choice], "!")
                else
                    println("Menu canceled.")
                end

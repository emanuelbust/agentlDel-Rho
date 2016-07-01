PYTHONS = $(wildcard *.py)
HELPS = $(wildcard *.help)
BACKUP = simBackup.tar.gz

$(BACKUP): $(PYTHONS) $(HELPS)
	tar -cvzf $(BACKUP) *.py *.help Makefile

.PHONY: clean
clean:
	rm -f $(BACKUP)
